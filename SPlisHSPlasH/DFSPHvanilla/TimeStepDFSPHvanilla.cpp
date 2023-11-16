#include "TimeStepDFSPHvanilla.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataDFSPHvanilla.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"


using namespace SPH;
using namespace std;
using namespace GenParam;


int TimeStepDFSPHvanilla::SOLVER_ITERATIONS_V = -1;
int TimeStepDFSPHvanilla::MAX_ITERATIONS_V = -1;
int TimeStepDFSPHvanilla::MAX_ERROR_V = -1;
int TimeStepDFSPHvanilla::USE_DIVERGENCE_SOLVER = -1;


TimeStepDFSPHvanilla::TimeStepDFSPHvanilla() :
	TimeStep(),
	m_simulationData()
{
	m_simulationData.init();

	m_counter = 0;
	m_iterationsV = 0;
	m_enableDivergenceSolver = true;
	m_maxIterationsV = 100;
	m_maxErrorV = static_cast<Real>(0.1);

	Simulation * sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		model->addField({ "factor", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getFactor(fluidModelIndex, i); } });
		model->addField({ "advected density", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getDensityAdv(fluidModelIndex, i); } });
		model->addField({ "p", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressure(fluidModelIndex, i); }, true });
		model->addField({ "p_v", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressure_V(fluidModelIndex, i); }, true });
		// model->addField({ "pressure acceleration", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureAccel(fluidModelIndex, i)[0]; } });
	}
}

TimeStepDFSPHvanilla::~TimeStepDFSPHvanilla(void)
{
	// remove all particle fields
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		model->removeFieldByName("factor");
		model->removeFieldByName("advected density");
		model->removeFieldByName("p");
		model->removeFieldByName("p_v");
		// model->removeFieldByName("pressure acceleration");
	}
}

void TimeStepDFSPHvanilla::initParameters()
{
	TimeStep::initParameters();

	 SOLVER_ITERATIONS_V = createNumericParameter("iterationsV", "Iterations (divergence)", &m_iterationsV);
	 setGroup(SOLVER_ITERATIONS_V, "Simulation|DFSPH");
	 setDescription(SOLVER_ITERATIONS_V, "Iterations required by the divergence solver.");
	 getParameter(SOLVER_ITERATIONS_V)->setReadOnly(true);

	 // MAX_ITERATIONS_V = createNumericParameter("maxIterationsV", "Max. iterations (divergence)", &m_maxIterationsV);
	 // setGroup(MAX_ITERATIONS_V, "Simulation|DFSPH");
	 // setDescription(MAX_ITERATIONS_V, "Maximal number of iterations of the divergence solver.");
	 // static_cast<NumericParameter<unsigned int>*>(getParameter(MAX_ITERATIONS_V))->setMinValue(1);
	 //
	 MAX_ERROR_V = createNumericParameter("maxErrorV", "Max. divergence error(%)", &m_maxErrorV);
	 setGroup(MAX_ERROR_V, "Simulation|DFSPH");
	 setDescription(MAX_ERROR_V, "Maximal divergence error (%).");
	 static_cast<RealParameter*>(getParameter(MAX_ERROR_V))->setMinValue(static_cast<Real>(1e-6));
	 //
	 // USE_DIVERGENCE_SOLVER = createBoolParameter("enableDivergenceSolver", "Enable divergence solver", &m_enableDivergenceSolver);
	 // setGroup(USE_DIVERGENCE_SOLVER, "Simulation|DFSPH");
	 // setDescription(USE_DIVERGENCE_SOLVER, "Turn divergence solver on/off.");
}

void TimeStepDFSPHvanilla::step()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	TimeManager* tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	performNeighborhoodSearch();

	// place for precomputing values to allow for avx-viscosity
#ifdef USE_PERFORMANCE_OPTIMIZATION
	//////////////////////////////////////////////////////////////////////////
	// precompute the values V_j * grad W_ij for all neighbors
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("precomputeValues")
	precomputeValues();
	STOP_TIMING_AVG
#endif

	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
		computeDensities(fluidModelIndex);
		computeDFSPHFactor(fluidModelIndex);
	}

	//////////////////////////////////////////////////////////////////////////
	// Divergence-Free Solver  -> make velocity field divergence-free
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("divergence-free solver");
	divergenceSolve();
	STOP_TIMING_AVG;

	//////////////////////////////////////////////////////////////////////////
	// update velocities using pressure accelerations to ensure that velocity field is divergence-free
	//////////////////////////////////////////////////////////////////////////
	// for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
	// 	FluidModel* fm = sim->getFluidModel(fluidModelIndex);
	// 	const unsigned int numParticles = fm->numActiveParticles();
	//
	// 	#pragma omp parallel default(shared)
	// 	{
	// 		#pragma omp for schedule(static)
	// 		for (int i = 0; i < numParticles; i++){
	// 			if (fm->getParticleState(i) == ParticleState::Active){
	// 				fm->getVelocity(i) += h * m_simulationData.getPressureAccel(fluidModelIndex, i);
	// 			}
	// 		}
	// 	}
	// }

	// reset the accelerations of the particles and add gravity
	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
		clearAccelerations(fluidModelIndex);
	}

	//////////////////////////////////////////////////////////////////////////
	// Non-pressure forces
	//////////////////////////////////////////////////////////////////////////
	// sim->computeNonPressureForces();
	// INFO: stored in acceleration array of fluid-model
	// note that neighbors and densities are already determined at this point
	// only viscosity (XSPH)
	// for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
	// 	for (int i = 0; i < Simulation::getCurrent()->getFluidModel(fluidModelIndex)->numActiveParticles(); i++){
	// 		computeViscosityForce(fluidModelIndex, i, h);
	// 	}
	// }

	//////////////////////////////////////////////////////////////////////////
	// advect velocities based on non-pressure forces incl. gravity
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = model->numActiveParticles();
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)
			for (int i = 0; i < numParticles; i++)
			{
				if (model->getParticleState(i) == ParticleState::Active)
				{
					Vector3r& vel = model->getVelocity(i);
					// accelerations already containing gravity and non-pressure accels
					vel += h * model->getAcceleration(i);
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// compute constant density source term
	//////////////////////////////////////////////////////////////////////////
	//for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
	//	FluidModel* fm = Simulation::getCurrent()->getFluidModel(fluidModelIndex);
	//	for (int i = 0; i < fm->numActiveParticles(); i++) {
	//        computeConstantDensitySourceTerm(fluidModelIndex, i, h);

	//		// Warm/Cold Start?
	//		Real& pressure = m_simulationData.getPressure(fluidModelIndex, i);
	//		pressure = 0.0;
	//    }
	//}


	//////////////////////////////////////////////////////////////////////////
	// Constant Density Solver
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("constant density solver");
	//pressureSolve();
	STOP_TIMING_AVG;


	//////////////////////////////////////////////////////////////////////////
	// update positions and velocities for pressure accelerations
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = model->numActiveParticles();

		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)
			for (unsigned int i = 0; i < numParticles; i++) {
				if (model->getParticleState(i) == ParticleState::Active)
				{
					Vector3r& xi = model->getPosition(i);
					Vector3r& vi = model->getVelocity(i);
					// vi += h * m_simulationData.getPressureAccel(fluidModelIndex, i);
					xi += h * vi;
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// emit new particles and perform an animation field step
	//////////////////////////////////////////////////////////////////////////
	sim->emitParticles();
	sim->animateParticles();

	//////////////////////////////////////////////////////////////////////////
	// Compute new time
	//////////////////////////////////////////////////////////////////////////
	tm->setTime(tm->getTime() + h);

	// update time step size
	sim->updateTimeStepSize();
}

void TimeStepDFSPHvanilla::computeDFSPHFactor(const unsigned int fluidModelIndex){
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	unsigned int numParticles = model->numActiveParticles();

	// TODO: OpenMP
	for(int i = 0; i < numParticles; i++){
		Real& factor = m_simulationData.getFactor(fluidModelIndex, i);

		const Real& rho_i = model->getDensity(i);
		const Real rho_i2 = rho_i * rho_i;
		const Real& density0_i = model->getDensity0();
		const Vector3r xi = model->getPosition(i);

		Real sumNormGradW = 0.0;
		Vector3r sumGrad = Vector3r::Zero();

		forall_fluid_neighbors(
			Vector3r mgw = fm_neighbor->getMass(neighborIndex) * sim->gradW(xi - xj);
			sumNormGradW += mgw.squaredNorm();
			sumGrad += mgw;
			);

		//////////////////////////////////////////////////////////////////////////
		// Boundary
		forall_boundary_neighbors(
			Vector3r mgw = bm_neighbor->getVolume(neighborIndex) * density0_i * sim->gradW(xi - xj);
			sumGrad += mgw;
			);
		//////////////////////////////////////////////////////////////////////////
		sumNormGradW += sumGrad.squaredNorm();

		factor = rho_i2 / (sumNormGradW + 0.0001);
	}
}

void TimeStepDFSPHvanilla::computeViscosityForce(const unsigned int fluidModelIndex, const unsigned i, const Real h) {
	Simulation* sim = Simulation::getCurrent();
	FluidModel* fm = Simulation::getCurrent()->getFluidModel(fluidModelIndex);
	Vector3r& accel_i = fm->getAcceleration(i);
	auto nFluids = Simulation::getCurrent()->numberOfFluidModels();
	const Vector3r &xi = fm->getPosition(i);

	Real factor = 0.01 / h;
	Vector3r sum = Vector3r::Zero();

	forall_fluid_neighbors(
		sum += (fm_neighbor->getMass(neighborIndex) / fm_neighbor->getDensity(neighborIndex)) * (fm_neighbor->getVelocity(neighborIndex) - fm->getVelocity(i)) * CubicKernel::W(xi - xj);
	);

	accel_i += factor * sum;
}

// void TimeStepDFSPHvanilla::pressureSolve(){
// 	Simulation* sim = Simulation::getCurrent();
// 	const unsigned int nFluids = sim->numberOfFluidModels();
// 	TimeManager* tm = TimeManager::getCurrent();
// 	const Real h = tm->getTimeStepSize();
//
// 	bool solved = false;
// 	m_iterations = 0;
//
// 	Real density_err = 0.0;
// 	Real avg_density_err = 0.0;
//
// 	while ((!solved || (m_iterations < m_minIterations)) && m_iterations < m_maxIterations){
// 		solved = true;
// 		avg_density_err = 0.0;
//
// 		//////////////////////////////////////////////////////////////////////////
// 		// compute pressure accelerations
// 		for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
// 			FluidModel* fm = sim->getFluidModel(fluidModelIndex);
// 			const unsigned int numParticles = fm->numActiveParticles();
//
// 			for (int i = 0; i < numParticles; i++) {
// 				computePressureAccel(fluidModelIndex, i, fm->getDensity0(), m_simulationData.getPressureRho2Data());
// 			}
// 		}
// 		//////////////////////////////////////////////////////////////////////////
//
// 		//////////////////////////////////////////////////////////////////////////
// 		// update pressure values
// 		for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
// 			FluidModel* fm = sim->getFluidModel(fluidModelIndex);
// 			const unsigned int numParticles = fm->numActiveParticles();
// 			const Real density0 = fm->getDensity0();
//
// 			density_err = 0.0;
// 			for(int i = 0; i < numParticles; i++) {
// 				const Real density_adv = m_simulationData.getDensityAdv(fluidModelIndex, i);
//
// 				const Real sourceTerm = m_simulationData.getSourceTerm(fluidModelIndex, i);
// 				Real& pressure_i = m_simulationData.getPressure(fluidModelIndex, i);
// 				Real& diag_i = m_simulationData.getDiagElement(fluidModelIndex, i);
// 				const Real aij_pj = compute_aij_pj(fluidModelIndex, i);
//
// 				// update pressure
// 				if (abs(diag_i) > m_eps){
// 					pressure_i += (static_cast<Real>(0.5) / diag_i) * (sourceTerm - aij_pj);
// 					// pressure_i += (static_cast<Real>(0.5) / diag_i) * ((density0 - density_adv) * (static_cast<Real>(1.0) / h)  - aij_pj);
// 				}
// 				else {
// 					pressure_i = 0.0;
// 				}
//
// 				// clamp to zero
// 				pressure_i = max(pressure_i, 0.0f);
// 				density_err -= min((sourceTerm - aij_pj), 0.0f) * h;
// 				// density_err -= min((density0 - density_adv) * (static_cast<Real>(1.0) / h) - aij_pj, 0.0f) * h;
// 			}
//
// 			avg_density_err = density_err / numParticles;
// 			solved = solved && (avg_density_err <= (m_maxErrorV * 0.01 * density0));
// 		 }
// 		//////////////////////////////////////////////////////////////////////////
// 		///
// 		m_iterations++;
// 	}
//
// 	// final acceleration update
// 	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
//		FluidModel* fm = sim->getFluidModel(fluidModelIndex);
//		const unsigned int numParticles = fm->numActiveParticles();
//
//		for (int i = 0; i < numParticles; i++) {
//			computePressureAccel(fluidModelIndex, i, fm->getDensity0(), m_simulationData.getPressureRho2Data());
//		}
//	}
// }

void TimeStepDFSPHvanilla::divergenceSolve()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();

	m_iterationsV = 0;
	bool solved = false;
	Real avg_divergence_error = 0.0;

	while((!solved || m_iterationsV < 1) && m_iterationsV < m_maxIterationsV) {
		solved = true;
		//////////////////////////////////////////////////////////////////////////
		// compute divergence source term and pressure values
		//////////////////////////////////////////////////////////////////////////
		for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
			Real divergence_error = 0.0;

			FluidModel* fm = sim->getFluidModel(fluidModelIndex);
			const unsigned int numParticles = fm->numActiveParticles();
			const Real density0 = fm->getDensity0();

			computeDivergenceSourceTerm(fluidModelIndex, h, divergence_error);
			computePressureV(fluidModelIndex, h);

			avg_divergence_error = divergence_error / numParticles;
			Real eta = m_maxErrorV; //* static_cast<Real>(0.01) * density0;
			solved = solved && (avg_divergence_error <= eta);
		}

		//////////////////////////////////////////////////////////////////////////
		// advect velocity
		//////////////////////////////////////////////////////////////////////////
		for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
			FluidModel* fm = sim->getFluidModel(fluidModelIndex);
			const Real density0 = fm->getDensity0();
			const unsigned int numParticles = fm->numActiveParticles();

			for (int i = 0; i < numParticles; i++) {
				Vector3r& vi = fm->getVelocity(i);
				const Vector3r ai_p = computePressureAccel(fluidModelIndex, i, density0, m_simulationData.getPressureRho2VData());
				vi += h * ai_p;
			}
		}

		m_iterationsV++;
	}

}

void TimeStepDFSPHvanilla::computeDivergenceSourceTerm(const unsigned int fluidModelIndex, const Real h, Real& divergence_error){
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int nFluids = sim->numberOfFluidModels();


	for(int i = 0; i < numParticles; i++){
		const Vector3r& xi = model->getPosition(i);

		Real& sourceTermDiv = m_simulationData.getSourceTermDiv(fluidModelIndex, i);
		Real denseDivVelocity = 0.0;

		forall_fluid_neighbors(
			denseDivVelocity += model->getMass(i) * (model->getVelocity(i) - fm_neighbor->getVelocity(neighborIndex)).dot(sim->gradW(xi - xj));
			);

		//////////////////////////////////////////////////////////////////////////
		// Boundary
		if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012) {
			forall_boundary_neighbors(
				denseDivVelocity += /* psi */ (bm_neighbor->getVolume(neighborIndex) * model->getDensity0()) * (model->getVelocity(i)).dot(sim->gradW(xi - xj));
				);
		}
		//////////////////////////////////////////////////////////////////////////

		sourceTermDiv = denseDivVelocity;
		divergence_error -= min(sourceTermDiv, 0.0f) * h;
	}
}

void TimeStepDFSPHvanilla::computePressureV(const unsigned int fluidModelIndex, const Real h){
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int numParticles = model->numActiveParticles();

	for (int i = 0; i < numParticles; i++){
		const Real sourceTerm = m_simulationData.getSourceTermDiv(fluidModelIndex, i);
		const Real factor = m_simulationData.getFactor(fluidModelIndex, i);
		Real& pressureV = m_simulationData.getPressure_V(fluidModelIndex, i);
		pressureV = 0.0;

		pressureV = (static_cast<Real>(1.0)/h) * factor * sourceTerm;
	}
}

void TimeStepDFSPHvanilla::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
	m_iterations = 0;
	m_iterationsV = 0;
}

void TimeStepDFSPHvanilla::performNeighborhoodSearch()
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

void TimeStepDFSPHvanilla::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepDFSPHvanilla::resize()
{
	m_simulationData.init();
}

/** Compute pressure accelerations using the current pressure values of the particles
 */
Vector3r TimeStepDFSPHvanilla::computePressureAccel(const unsigned int fluidModelIndex, const unsigned int i, const Real density0, const std::vector<std::vector<Real>>& pressureList, const bool applyBoundaryForces)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int nFluids = sim->numberOfFluidModels();

	Vector3r pressureAccel = Vector3r::Zero();

	const Real& p_i = pressureList[fluidModelIndex][i];
	const Real density_i = model->getDensity(i);
	const Real dpi = p_i / (density_i * density_i);

	const Vector3r& xi = model->getPosition(i);

	
	forall_fluid_neighbors(
		const Real& mass = fm_neighbor->getMass(neighborIndex);
		const Real& p_j = pressureList[pid][neighborIndex];
		const Real density_j = fm_neighbor->getDensity(neighborIndex);
		const Real dpj = p_j / (density_j * density_j);
		const Vector3r gradW = sim->gradW(xi - xj);

		pressureAccel -= mass * (dpi + dpj) * gradW;
	);

	//////////////////////////////////////////////////////////////////////////
	// Boundary Handling
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors(
			const Vector3r grad_p_j = /*psi*/ (bm_neighbor->getVolume(neighborIndex) * density0) * sim->gradW(xi - xj);
			const Vector3r a = (Real)1.0 * dpi * grad_p_j;
			pressureAccel -= a;
		);
	}
	//////////////////////////////////////////////////////////////////////////

	return pressureAccel;
}

void TimeStepDFSPHvanilla::compute_aii(const unsigned int fluidModelIndex, const unsigned int i, const Real h){
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int nFluids = sim->numberOfFluidModels();

	const Real& rho_i = model->getDensity(i);
	const Real& density0_i = model->getDensity0();
	const Vector3r xi = model->getPosition(i);

	Real sumNormGradW = 0.0;
	Vector3r sumGrad = Vector3r::Zero();

	forall_fluid_neighbors(
		Vector3r mgw = fm_neighbor->getMass(neighborIndex) * sim->gradW(xi - xj);
		sumNormGradW += mgw.squaredNorm();
		sumGrad += mgw;
	);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	forall_boundary_neighbors(
		Vector3r mgw = bm_neighbor->getVolume(neighborIndex) * density0_i * sim->gradW(xi - xj);
		sumGrad += mgw;
	);
	//////////////////////////////////////////////////////////////////////////


	sumNormGradW += sumGrad.squaredNorm();
	// m_simulationData.setDiagElement(fluidModelIndex, i, (-h / (rho_i * rho_i)) * sumNormGradW);
}

// Real TimeStepDFSPHvanilla::compute_aij_pj(const unsigned int fluidModelIndex, const unsigned int i)
// {
// 	Simulation* sim = Simulation::getCurrent();
// 	FluidModel* model = sim->getFluidModel(fluidModelIndex);
// 	const Real h = TimeManager::getCurrent()->getTimeStepSize();
// 	const unsigned int nFluids = sim->numberOfFluidModels();
// 	const Real density0 = model->getDensity0();
// 	//////////////////////////////////////////////////////////////////////////
// 	// Compute A*p which is the change of the density when applying the
// 	// pressure forces.
// 	//////////////////////////////////////////////////////////////////////////
// 	const Vector3r& xi = model->getPosition(i);
// 	const Vector3r& ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
// 	Real aij_pj = 0.0;
// 	//////////////////////////////////////////////////////////////////////////
// 	// Fluid
// 	//////////////////////////////////////////////////////////////////////////
// 	forall_fluid_neighbors(
// 		const Vector3r& aj = m_simulationData.getPressureAccel(pid, neighborIndex);
// 		aij_pj += fm_neighbor->getMass(neighborIndex) * (ai - aj).dot(sim->gradW(xi - xj));
// 		);
// 	//////////////////////////////////////////////////////////////////////////
// 	// Boundary
// 	//////////////////////////////////////////////////////////////////////////
// 	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
//	{
//		forall_boundary_neighbors(
//			aij_pj += bm_neighbor->getVolume(neighborIndex) * density0 * ai.dot(sim->gradW(xi - xj));
//		);
//	}
//
// 	return aij_pj * h;
// }


// void TimeStepDFSPHvanilla::computeConstantDensitySourceTerm (const unsigned int fluidModelIndex, const unsigned int i, const Real h){
//
// 	Simulation* sim = Simulation::getCurrent();
// 	const unsigned int nFluids = sim->numberOfFluidModels();
// 	FluidModel* model = sim->getFluidModel(fluidModelIndex);
// 	Real& sourceTerm = m_simulationData.getSourceTerm(fluidModelIndex, i);
// 	const Real density0 = model->getDensity0();
// 	const Real density = model->getDensity(i);
// 	const Vector3r& xi = model->getPosition(i);
// 	const Vector3r& vi = model->getVelocity(i);
//
// 	Real Drho_Dt = static_cast<Real>(0.0);
//
// 	forall_fluid_neighbors (
//         const Vector3r gradW = sim->gradW(xi - xj);
// 		Drho_Dt += model->getMass(i) * (vi - fm_neighbor->getVelocity(neighborIndex)).dot(gradW);
//     );
//
// 	//////////////////////////////////////////////////////////////////////////
// 	// Boundary
// 	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012){
// 		forall_boundary_neighbors(
// 			const Vector3r gradW = sim->gradW(xi - xj);
// 			Drho_Dt += (bm_neighbor->getVolume(neighborIndex) * density0) * vi.dot(gradW);
// 		);
// 	}
// 	//////////////////////////////////////////////////////////////////////////
//
// 	sourceTerm = 1.0 / h * (density0 - (density + h * Drho_Dt));
// }




