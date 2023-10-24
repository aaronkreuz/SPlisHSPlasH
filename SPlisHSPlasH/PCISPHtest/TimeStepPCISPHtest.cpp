#include "TimeStepPCISPHtest.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataPCISPHtest.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "../Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"

using namespace SPH;
using namespace std;

TimeStepPCISPHtest::TimeStepPCISPHtest() :
	TimeStep()
{
	m_simulationData.init();
	m_counter = 0;
	m_minIterations = 3;
	m_maxIterations = 100;

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		
		// INFO: addField (FieldDescription::name, FieldDescription::type, FieldDescription::getFunction)
		// INFO: FieldDescription::getFunction returns a pointer to the function that returns the value of the field
	
		model->addField({ "density_adv", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getDensityAdv(fluidModelIndex, i); } });
		model->addField({ "pressure", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressure(fluidModelIndex, i); } });
		model->addField({ "pressure acceleration", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureAccel(fluidModelIndex, i)[0]; } });
	}
}

TimeStepPCISPHtest::~TimeStepPCISPHtest(void)
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);

		model->removeFieldByName("density_adv");
		model->removeFieldByName("pressure");
		model->removeFieldByName("pressure acceleration");
	}
}

void TimeStepPCISPHtest::step(){
	const auto nFluids = Simulation::getCurrent()->numberOfFluidModels();
	Simulation* sim = Simulation::getCurrent();
	// TimeManager* tm = TimeManager::getCurrent();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();

	// clear accelerations already adds gravity
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
		clearAccelerations(fluidModelIndex);

	performNeighborhoodSearch();

#ifdef USE_PERFORMANCE_OPTIMIZATION
	precomputeValues();
#endif

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		computeVolumeAndBoundaryX();
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		computeDensityAndGradient();

	for (auto fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
		computeDensities(fluidModelIndex);
	}
	//////////////////////////////////////////////////////////////////////////

	sim->computeNonPressureForces();

	sim->updateTimeStepSize();

	// TODO: Multi-Threading
	for (auto fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
		FluidModel* fm = sim->getFluidModel(fluidModelIndex);
		for(auto i = 0; i < fm->numActiveParticles(); i++){

			Vector3r& vel = fm->getVelocity(i);
			const Vector3r& accel = fm->getAcceleration(i);
			if (fm->getParticleState(i) == ParticleState::Active)
				vel += h * accel;

			m_simulationData.getPressure(fluidModelIndex, i) = 0.0;
			m_simulationData.getPressureAccel(fluidModelIndex, i).setZero();
		}
	}

	solvePressure();

	timeIntegration();
}


void TimeStepPCISPHtest::solvePressure(){
	Simulation* sim = Simulation::getCurrent();
	auto nFluids = Simulation::getCurrent()->numberOfFluidModels();
	int solverIterations = 0;
	bool solved = false;

	Real avg_density_err = 0.0;

	////////////////////////////////////////////////////////////////////
	// pressure solver loop
	///////////////////////////////////////////////////////////////////
	while (solverIterations < m_minIterations || (!solved && (solverIterations < m_maxIterations))){

		solved = true;

#pragma omp parallel default(shared)
		{
#pragma omp for schedule(static)
			for(int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
				FluidModel* fm = Simulation::getCurrent()->getFluidModel(fluidModelIndex);
				for(auto i = 0; i < fm->numActiveParticles(); i++){
					Vector3r& predV = m_simulationData.getPredictedVelocity(fluidModelIndex, i);
					Vector3r& predX = m_simulationData.getPredictedPosition(fluidModelIndex, i);
					const Vector3r& vel = fm->getVelocity(i);
					const Vector3r& pos = fm->getPosition(i);

					predV = vel + TimeManager::getCurrent()->getTimeStepSize() * (fm->getAcceleration(i) + m_simulationData.getPressureAccel(fluidModelIndex, i));
					predX = pos + TimeManager::getCurrent()->getTimeStepSize() * predV;
					if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
						computeVolumeAndBoundaryX(fluidModelIndex, i, predX);
					else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
						computeDensityAndGradient(fluidModelIndex, i, predX);
				}
			}
		}

		for(unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
			const Real density0 = sim->getFluidModel(fluidModelIndex)->getDensity0();
			avg_density_err = 0.0;

			pressureSolverIteration(fluidModelIndex, avg_density_err);


			// Maximal allowed density fluctuation
			const Real eta = m_maxError * static_cast<Real>(0.01) * density0;  // maxError is given in percent
			solved = solved && (avg_density_err <= eta);
		}

		solverIterations++;
	}
}

void TimeStepPCISPHtest::pressureSolverIteration(int fluidModelIndex, Real& avg_density_err){
	Simulation* sim = Simulation::getCurrent();
	FluidModel* fm = sim->getFluidModel(fluidModelIndex);
	int numParticles = fm->numActiveParticles();
	auto nFluids = sim->numberOfFluidModels();
	auto nBoundaries = sim->numberOfBoundaryModels();
	const auto density0 = fm->getDensity0();

	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real h2 = h*h;
	const Real invH2 = static_cast<Real>(1.0) / h2;

	Real sumError = 0.0;


// Predict density
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)
		for (int i = 0; i < numParticles; i++)
		{
			const Vector3r &pred_xi = m_simulationData.getPredictedPosition(fluidModelIndex, i);
			Real &densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
			densityAdv = fm->getVolume(i) * sim->W_zero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Vector3r & pred_xj = m_simulationData.getPredictedPosition(pid, neighborIndex);
				densityAdv += fm_neighbor->getVolume(neighborIndex) * sim->W(pred_xi - pred_xj);
				);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					densityAdv += bm_neighbor->getVolume(neighborIndex) * sim->W(pred_xi - xj);
					);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					densityAdv += rho;
					);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					densityAdv += Vj * sim->W(pred_xi - xj);
					);
			}

			//densityAdv = max(densityAdv, static_cast<Real>(1.0));
			const Real density_err = density0 * (densityAdv - static_cast<Real>(1.0));
			Real &pressure = m_simulationData.getPressure(fluidModelIndex, i);
			pressure += invH2 * m_simulationData.getPcisphFactor(fluidModelIndex) * (densityAdv - static_cast<Real>(1.0));
			pressure = max(pressure, static_cast<Real>(0.0));

			#pragma omp atomic
			sumError += density_err;
		}
	}

	avg_density_err = sumError / numParticles;

	// for(auto i = 0; i < numParticles; i++){
	// 	//////////////////////////////////////////////////////////////////////////
	// 	// compute advected density
	// 	//////////////////////////////////////////////////////////////////////////
	// 	const Vector3r& pred_xi = m_simulationData.getPredictedPosition(fluidModelIndex, i);
	// 	auto& advDensity = m_simulationData.getDensityAdv(fluidModelIndex, i);

	// 	advDensity = fm->getVolume(i) * sim->W_zero();
	// 	forall_fluid_neighbors(
	// 		const Vector3r& pred_xj = m_simulationData.getPredictedPosition(pid, neighborIndex);
	// 		advDensity += fm_neighbor->getVolume(neighborIndex) * sim->W(pred_xi - pred_xj);
	// 	);

	// 	//////////////////////////////////////////////////////////////////////////
	// 	// Boundary
	// 	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	// 	{
	// 		forall_boundary_neighbors(
	// 			advDensity += bm_neighbor->getVolume(neighborIndex) * sim->W(pred_xi - xj);
	// 			);
	// 	}
	// 	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	// 	{
	// 		forall_density_maps(
	// 			advDensity += rho;
	// 			);
	// 	}
	// 	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	// 	{
	// 		forall_volume_maps(
	// 			advDensity += Vj * sim->W(pred_xi - xj);
	// 			);
	// 	}
	// 	//////////////////////////////////////////////////////////////////////////

	// 	//////////////////////////////////////////////////////////////////////////
	// 	// predict density variation
	// 	//////////////////////////////////////////////////////////////////////////
	// 	// solve for Volume (multiply with rest density)
	// 	const Real densityVariation = (advDensity * fm->getDensity0()) - fm->getDensity0(); // i.e. error
	// 	sumError += densityVariation;

	// 	//////////////////////////////////////////////////////////////////////////
	// 	// update pressure
	// 	//////////////////////////////////////////////////////////////////////////
	// 	auto& pressure = m_simulationData.getPressure(fluidModelIndex, i);
	// 	const Real pcisphFactor = m_simulationData.getPcisphFactor(fluidModelIndex);
	// 	pressure += invH2 * pcisphFactor * densityVariation;
	// 	pressure = max(pressure, static_cast<Real>(0.0));
	// }

		// for( int i = 0; i < numParticles; i++){
		// //////////////////////////////////////////////////////////////////////////
		// // compute pressure acceleration
		// //////////////////////////////////////////////////////////////////////////
		// Vector3r &pred_xi = m_simulationData.getPredictedPosition(fluidModelIndex, i);
		// const Real &pressure = m_simulationData.getPressure(fluidModelIndex, i);
		// Vector3r pressureAccel = m_simulationData.getPressureAccel(fluidModelIndex, i);
		// pressureAccel.setZero();

		// const Real& m_i = fm->getMass(i);
		// const Real m2 = m_i * m_i;
		// const Real rho2 = fm->getDensity0() * fm->getDensity0();
		// Vector3r sumGradKernel = Vector3r::Zero();
		// const Vector3r& xi = fm->getPosition(i);

		// forall_fluid_neighbors(
		// 	const Vector3r& pred_xj = m_simulationData.getPredictedPosition(pid, neighborIndex);
		// 	sumGradKernel += sim->gradW(pred_xi - pred_xj);
		// );

		// pressureAccel = static_cast<Real>(2.0) * m2 *( pressure / rho2) * sumGradKernel;

		// //////////////////////////////////////////////////////////////////////////
		// // Boundary
		// if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		// {
		// 	forall_boundary_neighbors(
		// 		const Vector3r a = bm_neighbor->getVolume(neighborIndex) * pressure * sim->gradW(xi - xj);
		// 		pressureAccel -= a;
		// 		bm_neighbor->addForce(xj, fm->getMass(i) * a);
		// 		);
		// }
		// else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		// {
		// 	forall_density_maps(
		// 		const Vector3r a = -pressure * gradRho;
		// 		pressureAccel -= a;
		// 		bm_neighbor->addForce(xj, fm->getMass(i) * a);
		// 		);
		// }
		// else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		// {
		// 	forall_volume_maps(
		// 		const Vector3r a = Vj * pressure * sim->gradW(xi - xj);
		// 		pressureAccel -= a;
		// 		bm_neighbor->addForce(xj, fm->getMass(i) * a);
		// 		);
		// }
		// //////////////////////////////////////////////////////////////////////////
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)
			for (int i = 0; i < (int)numParticles; i++)
			{
				const Vector3r &xi = fm->getPosition(i);

				Vector3r &ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
				ai.setZero();

				const Real dpi = m_simulationData.getPressure(fluidModelIndex, i);

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				forall_fluid_neighbors(
					// Pressure
					const Real dpj = m_simulationData.getPressure(pid, neighborIndex);
					ai -= fm_neighbor->getVolume(neighborIndex) * (dpi + (fm_neighbor->getDensity0() / density0) * dpj) * sim->gradW(xi - xj);
					);

				//////////////////////////////////////////////////////////////////////////
				// Boundary
				//////////////////////////////////////////////////////////////////////////
				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
				{
					forall_boundary_neighbors(
						// Pressure
						const Vector3r a = bm_neighbor->getVolume(neighborIndex) * dpi * sim->gradW(xi - xj);
						ai -= a;
						bm_neighbor->addForce(xj, fm->getMass(i) * a);
						);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
				{
					forall_density_maps(
						const Vector3r a = -dpi * gradRho;
						ai -= a;
						bm_neighbor->addForce(xj, fm->getMass(i) * a);
						);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
				{
					forall_volume_maps(
						const Vector3r a = Vj *dpi * sim->gradW(xi - xj);
						ai -= a;
						bm_neighbor->addForce(xj, fm->getMass(i) * a);
						);
				}
			}
		}
}

void TimeStepPCISPHtest::timeIntegration(){
	auto* sim = Simulation::getCurrent();
	auto nFluids = sim->numberOfFluidModels();
	auto* tm = TimeManager::getCurrent();
	Real h = tm->getTimeStepSize();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const int numParticles = (int)model->numActiveParticles();

#pragma omp parallel default(shared)
		{
#pragma omp for schedule(static)
			for (int i = 0; i < numParticles; i++)
			{
				if (model->getParticleState(i) == ParticleState::Active)
				{
					const Vector3r &accel = m_simulationData.getPressureAccel(fluidModelIndex, i) + model->getAcceleration(i);
					Vector3r &x = model->getPosition(i);
					Vector3r &v = model->getVelocity(i);
					v += h * accel;
					x += h * v;
				}
			}
		}
	}

	sim->emitParticles();
	sim->animateParticles();

	// Compute new time
	tm->setTime(tm->getTime() + h);
}

void TimeStepPCISPHtest::step2()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	TimeManager* tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	// clear accelerations already adds gravity
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
		clearAccelerations(fluidModelIndex);

	// 1) compute current neighborhood for particlesnumParticles
	performNeighborhoodSearch();

#ifdef USE_PERFORMANCE_OPTIMIZATION
	precomputeValues();
#endif

	// Required for the non-pressure force computations
	// TOOO: Understand this
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		computeVolumeAndBoundaryX();
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		computeDensityAndGradient();

	// 2) reconstruct densities at particle positions for every fluid model
	for (auto i = 0; i < nFluids; i++) {
		// densitites are stored in corresponding fluid model
		computeDensities(i);
	}

	// DEBUG
	// LOG_INFO << "densities computed";

	////////////////////////////////////////////////////////////////////////////
	// 3) COMPUTE NON-PRESSURE FORCES
	////////////////////////////////////////////////////////////////////////////
	// stored in the acceleration field of the corresponding fluid model
	sim->computeNonPressureForces();

	// DEBUG
	// LOG_INFO << "non-pressure forces computed";

	// TODO: 3.1) already update velocity with nonPressureForce accelerations ???

	// 3.2) pressure solver loop initialization
	for(auto fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
		FluidModel* fModel = sim->getFluidModel(fluidModelIndex);

		// INFO: only iterate over active particles
		const int numParticles = (int)fModel->numActiveParticles();
		for(auto i = 0; i < numParticles; i++) {
			Vector3r& vel = fModel->getVelocity(i);
			const Vector3r& accel = fModel->getAcceleration(fluidModelIndex);

			// TODO: add non-pressure accelerations already here? Should they then be omitted in the prediction-correction loop?
			if (fModel->getParticleState(fluidModelIndex) == ParticleState::Active)
				vel += h * accel;

			// set pressure and pressure acceleration to 0
			m_simulationData.getPressure(fluidModelIndex, i) = 0.0;
			m_simulationData.getPressureAccel(fluidModelIndex, i).setZero();
		}
	}

	// DEBUG
	// LOG_INFO << "pressure and pressure acceleration set to 0";

	////////////////////////////////////////////////////////////////////
	// 4) PRESSURE SOLVER - PREDICTION-CORRECTION SCHEME
	///////////////////////////////////////////////////////////////////
	m_iterations = 0;

	bool pressureSolved = false;
	LOG_INFO << "starting pressure solver loop";

	while (m_iterations < m_maxIterations && ((m_iterations < m_minIterations) || !pressureSolved)) {
		const unsigned int nBoundaries = sim->numberOfBoundaryModels();

		// 4.1) predict velocities and positions
		for (auto fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
			FluidModel* fModel = sim->getFluidModel(fluidModelIndex);
			// only active particles
			const int numParticles = (int)fModel->numActiveParticles();


			// define parallel region
			#pragma omp parallel default(shared)
			{
				// loop construct and schedule clause (static -> chunks are assigned to threads statically in Round-Robin manner in order of thread number)
				#pragma omp for schedule(static) nowait
				for (auto i = 0; i < numParticles; i++) {

					if (fModel->getParticleState(i) != ParticleState::Active)
						continue;

					auto& predV = m_simulationData.getPredictedVelocity(fluidModelIndex, i);
					auto& predX = m_simulationData.getPredictedPosition(fluidModelIndex, i);
					const auto& velocity = fModel->getVelocity(i);
					const auto& position = fModel->getPosition(i);

					predV = velocity + h * (fModel->getAcceleration(i) + m_simulationData.getPressureAccel(fluidModelIndex, i));
					predX = position + h * predV;
					if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
						computeVolumeAndBoundaryX(fluidModelIndex, i, predX);
					else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
						computeDensityAndGradient(fluidModelIndex, i, predX);
				}
			}
        }

		// DEBUG
		// LOG_INFO << "Pressure Solver | " << "predicted velocities and positions";

		// 4.2) Recompute smoothing kernel values
		// -> not necessary!

		/////////////////////////////////////////////////////////////////////
        // 4.3) pressure solver iteration
		////////////////////////////////////////////////////////////////////
		for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
			FluidModel* fm = sim->getFluidModel(fluidModelIndex);
			// only active particles
			const int numParticles = (int)fm->numActiveParticles();

			Real error_sum = 0.0;

			for (auto i = 0; i < numParticles; i++) {
				// 4.3a) Compute predicted density
				// INFO: Predicted (advected) density is stored in the simulation data
				auto& advDensity = m_simulationData.getDensityAdv(fluidModelIndex, i);
				auto& pred_xi = m_simulationData.getPredictedPosition(fluidModelIndex, i);
				const Vector3r &xi = fm->getPosition(i);

				// INFO: Naming of surrounding indices is fixed. I.e. fluidModelIndex, i, nFluids...
				// INFO: pid is the index of the fluid model of neighor particle, neighborIndex is the index of the neighbor particle
				// recomputing density at t+1 with standard SPH approximation
				forall_fluid_neighbors(
					const Vector3r& pred_xj = m_simulationData.getPredictedPosition(pid, neighborIndex);
					advDensity += fm_neighbor->getMass(neighborIndex) * sim->W(pred_xi - pred_xj);
				);
				//////////////////////////////////////////////////////////////////////////
				// Boundary handling
				//////////////////////////////////////////////////////////////////////////
				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
				{
					forall_boundary_neighbors(
						advDensity += bm_neighbor->getVolume(neighborIndex) * sim->W(pred_xi - xj);
						);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
				{
					forall_density_maps(
						advDensity += rho;
						);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
				{
					forall_volume_maps(
						advDensity += Vj * sim->W(pred_xi - xj);
						);
				}

				// advected density error
				const auto density_error = fm->getDensity0() - advDensity;

				error_sum += density_error;

				// 4.3b) Compute pressure
				// INFO: Pressure is stored in the simulation data
				auto& pressure = m_simulationData.getPressure(fluidModelIndex, i);
				pressure += m_simulationData.getPcisphFactor(fluidModelIndex) * (fm->getDensity0() - advDensity);

				// 4.3c) Compute pressure acceleration (force)
				// INFO: Pressure acceleration is stored in the simulation data
				auto& pressureAccel = m_simulationData.getPressureAccel(fluidModelIndex, i);
				pressureAccel.setZero();

				// equation from Solenthaler and Pajarola [SP09] PCISPH
				const Real pressureAccelFactor = (fm->getMass(i) * static_cast<Real>(2.0) * pressure) / (fm->getDensity0() * fm->getDensity0());
				Vector3r sumGradW = Vector3r::Zero();

				forall_fluid_neighbors(
					sumGradW += sim->gradW(m_simulationData.getPredictedPosition(pid, neighborIndex) - m_simulationData.getPredictedPosition(fluidModelIndex, i));
				);
				// LOG_INFO << "sumGradW computed for particle " << i << "|" << numParticles;


				// LOG_INFO << "pressureAccelFactor: " << pressureAccelFactor;
				// LOG_INFO << "sumGradW: " << sumGradW;
				pressureAccel = -pressureAccelFactor * sumGradW;
				// LOG_INFO << "pressureAccel: " << pressureAccel;

				//////////////////////////////////////////////////////////////////////////
				// Boundary
				//////////////////////////////////////////////////////////////////////////
				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
				{
					forall_boundary_neighbors(
						// Pressure
						const Vector3r a = bm_neighbor->getVolume(neighborIndex) * pressure * sim->gradW(xi - xj);
						pressureAccel -= a;
						bm_neighbor->addForce(xj, fm->getMass(i) * a);
						);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
				{
					forall_density_maps(
						const Vector3r a = -pressure * gradRho;
						pressureAccel -= a;
						bm_neighbor->addForce(xj, fm->getMass(i) * a);
						);
				}
				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
				{
					forall_volume_maps(
						const Vector3r a = Vj *pressure * sim->gradW(xi - xj);
						pressureAccel -= a;
						bm_neighbor->addForce(xj, fm->getMass(i) * a);
						);
				}
			}

			// LOG_INFO << "computed pressure, pressure acceleration and density for all particles of fluid " << fluidModelIndex;

			// 4.3c) Average error (per fluid model)
			// INFO: error_sum is the sum of all errors of all particles of all fluid models
			// INFO: error_sum is divided by the number of particles of all fluid models
			Real avg_error = error_sum / fm->numParticles();

			// Condition from Ihmsen&Tescher 2010
			pressureSolved = pressureSolved && (abs(avg_error) < (fm->getDensity0() * static_cast<Real>(0.02)) );

			m_iterations++;
		}


	}

	// 5) Update velocities and positions
	for (auto fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
		FluidModel* fm = sim->getFluidModel(fluidModelIndex);
		for (auto i = 0; i < fm->numParticles(); i++) {

			auto& position = fm->getPosition(i);
			auto& velocity = fm->getVelocity(i);

			// TODO: Add again the non-pressure acceleration?
			// fm->getAcceleration() contains the non-pressure accelerations
			velocity += h * (m_simulationData.getPressureAccel(fluidModelIndex, i) + fm->getAcceleration(i));
			position += h * velocity;
		}
	}
	
}

void TimeStepPCISPHtest::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
}

void TimeStepPCISPHtest::performNeighborhoodSearch()
{
	if (Simulation::getCurrent()->zSortEnabled())
	{
		if (m_counter % 500 == 0)
		{
			Simulation::getCurrent()->performNeighborhoodSearchSort();
			m_simulationData.performNeighborhoodSearchSort();
		}
		m_counter++;
	}

	Simulation::getCurrent()->performNeighborhoodSearch();
}

void TimeStepPCISPHtest::emittedParticles(FluidModel* model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepPCISPHtest::resize()
{
	m_simulationData.init();
}
