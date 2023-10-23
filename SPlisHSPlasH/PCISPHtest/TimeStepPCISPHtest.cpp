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
	std::cout << "init simulation data done" << std::endl;
	m_counter = 0;
	m_minIterations = 3;

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

void TimeStepPCISPHtest::step()
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
	// LOG_INFO << "starting pressure solver loop";

	while (m_iterations < m_maxIterations && ((m_iterations < m_minIterations) || !pressureSolved)) {
		const unsigned int nBoundaries = sim->numberOfBoundaryModels();

		// 4.1) predict velocities and positions
		for (auto fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
			FluidModel* fModel = sim->getFluidModel(fluidModelIndex);
			// only active particles
			const int numParticles = (int)fModel->numActiveParticles();

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
				LOG_INFO << "sumGradW computed for particle " << i << "|" << numParticles;


				// LOG_INFO << "pressureAccelFactor: " << pressureAccelFactor;
				// LOG_INFO << "sumGradW: " << sumGradW;
				pressureAccel = -pressureAccelFactor * sumGradW;
				LOG_INFO << "pressureAccel: " << pressureAccel;

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

			LOG_INFO << "computed pressure, pressure acceleration and density for all particles of fluid " << fluidModelIndex;
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
