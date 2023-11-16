#include "TimeStepDFSPHbubble.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
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


int TimeStepDFSPHbubble::SOLVER_ITERATIONS_V = -1;
int TimeStepDFSPHbubble::MAX_ITERATIONS_V = -1;
int TimeStepDFSPHbubble::MAX_ERROR_V = -1;
int TimeStepDFSPHbubble::USE_DIVERGENCE_SOLVER = -1;
int TimeStepDFSPHbubble::USE_COHESION_FORCE = -1;
int TimeStepDFSPHbubble::USE_DRAG_FORCE_ON_LIQ = -1;
int TimeStepDFSPHbubble::USE_DRAG_FORCE_ON_AIR = -1;
int TimeStepDFSPHbubble::USE_VISCOSITY = -1;
int TimeStepDFSPHbubble::USE_BOUYANCY = -1;
int TimeStepDFSPHbubble::USE_SURFACE_TENSION = -1;
int TimeStepDFSPHbubble::MAX_K_BOUYANCY = -1;
int TimeStepDFSPHbubble::MIN_BOUYANCY = -1;

int TimeStepDFSPHbubble::COHESION_FORCE_TYPE = -1;
int TimeStepDFSPHbubble::ENUM_COHESION_FORCE_IHMSEN = -1;
int TimeStepDFSPHbubble::ENUM_COHESION_FORCE_SURFACE_TENSION = -1;


TimeStepDFSPHbubble::TimeStepDFSPHbubble() :
	TimeStep(),
	m_simulationData()
{
	m_simulationData.init();

	m_counter = 0;
	m_iterationsV = 0;

	m_enableDivergenceSolver = true;
	m_enableCohesionForce = true;
	m_enableDragForceOnAir = true;
	m_enableDragForceOnLiq = true;
	m_enableViscosity = true;
	m_enableBouyancy = true;
	m_enableSurfaceTension = true;

	m_minBouyancy = static_cast<Real>(1.40);

	m_cohesionForce = 0; //. i.e. standard cohesion force described by Ihmsen et al. in "Animation of air bubbles with SPH"

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
		model->addField({ "pressure acceleration", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureAccel(fluidModelIndex, i)[0]; } });
	}
}

TimeStepDFSPHbubble::~TimeStepDFSPHbubble(void)
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
		model->removeFieldByName("pressure acceleration");
	}
}

void TimeStepDFSPHbubble::initParameters()
{
	TimeStep::initParameters();

	 USE_COHESION_FORCE = createBoolParameter("enableCohesionForce", "Enable cohesion force", &m_enableCohesionForce);
	 setGroup(USE_COHESION_FORCE, "Simulation|Forces");
	 setDescription(USE_COHESION_FORCE, "turn cohesion force on/off.");

	 USE_DRAG_FORCE_ON_AIR = createBoolParameter("enableDragForceAir", "Enable drag force on air", &m_enableDragForceOnAir);
	 setGroup(USE_DRAG_FORCE_ON_AIR, "Simulation|Forces");
	 setDescription(USE_DRAG_FORCE_ON_AIR, "turn drag force on air on/off.");

	 USE_DRAG_FORCE_ON_LIQ = createBoolParameter("enableDragForceLiquid", "Enable drag force on liquid", &m_enableDragForceOnLiq);
	 setGroup(USE_DRAG_FORCE_ON_LIQ, "Simulation|Forces");
	 setDescription(USE_DRAG_FORCE_ON_LIQ, "turn drag force on liquid on/off.");

	 USE_VISCOSITY = createBoolParameter("enableViscosity", "Enable viscosity", &m_enableViscosity);
	 setGroup(USE_VISCOSITY, "Simulation|Forces");
	 setDescription(USE_VISCOSITY, "turn viscosity on/off.");

	 USE_BOUYANCY = createBoolParameter("enableBouyancy", "Enable bouyancy", &m_enableBouyancy);
	 setGroup(USE_BOUYANCY, "Simulation|Forces");
	 setDescription(USE_BOUYANCY, "turn bouyancy on/off.");

	 USE_SURFACE_TENSION = createBoolParameter("enableSurfaceTension", "Enable surface tension", &m_enableSurfaceTension);
	 setGroup(USE_SURFACE_TENSION, "Simulation|Forces");
	 setDescription(USE_SURFACE_TENSION, "turn surface tension on/off.");

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

	 USE_DIVERGENCE_SOLVER = createBoolParameter("enableDivergenceSolver", "Enable divergence solver", &m_enableDivergenceSolver);
	 setGroup(USE_DIVERGENCE_SOLVER, "Simulation|DFSPH");
	 setDescription(USE_DIVERGENCE_SOLVER, "Turn divergence solver on/off.");

	 MIN_BOUYANCY = createNumericParameter("minBouyancy", "Min. bouyancy", &m_minBouyancy);
	 setGroup(MIN_BOUYANCY, "Simulation|BUBBLE");
	 setDescription(MIN_BOUYANCY, "Minimal bouyancy coefficient.");
	 static_cast<RealParameter*>(getParameter(MIN_BOUYANCY))->setMinValue(static_cast<Real>(1.0));

	 MAX_K_BOUYANCY = createNumericParameter("max_KBouyancy", "Max. K Bouyancy", &m_kmax);
	 setGroup(MAX_K_BOUYANCY, "Simulation|BUBBLE");
	 setDescription(MAX_K_BOUYANCY, "Maximal k bouyancy coefficient.");
	 static_cast<RealParameter*>(getParameter(MAX_K_BOUYANCY))->setMinValue(static_cast<Real>(1.0));

	 COHESION_FORCE_TYPE = createEnumParameter("cohesionForceType", "Cohesion force", &m_cohesionForce);
	 setGroup(COHESION_FORCE_TYPE, "Simulation|BUBBLE");
	 setDescription(COHESION_FORCE_TYPE, "Method for the cohesion force computation.");
	 EnumParameter *enumParam = static_cast<EnumParameter*>(getParameter(COHESION_FORCE_TYPE));
	 enumParam->addEnumValue("Ihmsen", ENUM_COHESION_FORCE_IHMSEN);
	 enumParam->addEnumValue("SurfaceTension", ENUM_COHESION_FORCE_SURFACE_TENSION);
}

void TimeStepDFSPHbubble::step()
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

	// compute densities separately for liquid and air
	// they don't see each other regarding densities -> coupling only via drag force!
	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
		FluidModel* fm = sim->getFluidModel(fluidModelIndex);

		// WARNING: JUST FOR TESTING
		// if(fm->getId() == "Liquid"){
		// 	computeDensitiesForLiquidPhase(fluidModelIndex);
		// }
		// else { // Air
		// 	computeDensitiesForSamePhase(fluidModelIndex);
		// }

		computeDensitiesForSamePhase(fluidModelIndex);
	}

	//////////////////////////////////////////////////////////////////////////
	// compute diagonal matrix elements -> same phase only (+ boundary)
	//////////////////////////////////////////////////////////////////////////
	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
		#pragma omp parallel default(shared)
		{
			FluidModel* fm = sim->getFluidModel(fluidModelIndex);
			const unsigned int numParticles = fm->numActiveParticles();

			#pragma for schedule(static)
			for (int i = 0; i < numParticles; i++) {
				compute_aii(fluidModelIndex, i, h);
			}
		}
	}

	if(m_enableDivergenceSolver){
		//////////////////////////////////////////////////////////////////////////
		// compute divergence source term
		//////////////////////////////////////////////////////////////////////////
		 for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
			FluidModel* fm = sim->getFluidModel(fluidModelIndex);
			#pragma omp parallel default(shared)
			{
				#pragma omp for schedule(static)
				for (int i = 0; i < fm->numActiveParticles(); i++) {
					computeDivergenceSourceTerm(fluidModelIndex, i, h);

					// Warm/Cold Start?
					Real& pressureV = m_simulationData.getPressure_V(fluidModelIndex, i);
					pressureV = 0.0;
				}
			}
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
		for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
			FluidModel* fm = sim->getFluidModel(fluidModelIndex);
			const unsigned int numParticles = fm->numActiveParticles();

			#pragma omp parallel default(shared)
			{
				#pragma omp for schedule(static)
				for (int i = 0; i < numParticles; i++){
					if (fm->getParticleState(i) == ParticleState::Active){
						fm->getVelocity(i) += h * m_simulationData.getPressureAccel(fluidModelIndex, i);
					}
				}
			}
		}
	}

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

	//////////////////////////////////////////////////////////////////////////
	// Viscosity XSPH
	//////////////////////////////////////////////////////////////////////////
	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
		FluidModel* fm = sim->getFluidModel(fluidModelIndex);

		// skip air -> viscosity not applied to air
		if(fm->getId() == "Air"){
			continue;
		}

		for (int i = 0; i < Simulation::getCurrent()->getFluidModel(fluidModelIndex)->numActiveParticles(); i++){
			computeViscosityForce(fluidModelIndex, i, h);
		}
	 }

	 //////////////////////////////////////////////////////////////////////////
	 // Surface Tension Force
	 //////////////////////////////////////////////////////////////////////////
	 if (m_enableSurfaceTension){
		computeSurfaceTensionForce(1, h);
	 }

	 //////////////////////////////////////////////////////////////////////////
	 // Drag Forces -> Two Way Coupling
	 //////////////////////////////////////////////////////////////////////////
	 // Force acting on Air particles
	 if (m_enableDragForceOnAir)
		computeDragForce(0, h);
	 // Force acting on Liquid particles
	 if (m_enableDragForceOnLiq)
		computeDragForce(1, h);


	 //////////////////////////////////////////////////////////////////////////
	 // Buoyancy Force -> only acting on Air particles
	 //////////////////////////////////////////////////////////////////////////
	 if (m_enableBouyancy)
		computeBouyancyForce(0, h);

	 //////////////////////////////////////////////////////////////////////////
	 // Cohesion Force -> only actiong on Air particles
	 //////////////////////////////////////////////////////////////////////////
	 if (m_enableCohesionForce){
		if(m_cohesionForce == ENUM_COHESION_FORCE_IHMSEN){
			computeCohesionForce(0,h);
		}
		// else if ..
	 }

	//////////////////////////////////////////////////////////////////////////
	// adapt time-step size
	//////////////////////////////////////////////////////////////////////////
	 sim->updateTimeStepSize();

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
	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
		FluidModel* fm = Simulation::getCurrent()->getFluidModel(fluidModelIndex);
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)
			for (int i = 0; i < fm->numActiveParticles(); i++) {
				computeConstantDensitySourceTerm(fluidModelIndex, i, h);

				// Warm/Cold Start?
				Real& pressure = m_simulationData.getPressure(fluidModelIndex, i);
				pressure = 0.0;
			}
		}
	}


	//////////////////////////////////////////////////////////////////////////
	// Constant Density Solver
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("constant density solver");
	pressureSolve();
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
					vi += h * m_simulationData.getPressureAccel(fluidModelIndex, i);
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
}

void TimeStepDFSPHbubble::pressureSolve(){
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	TimeManager* tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	bool solved = false;
	m_iterations = 0;

	Real density_err = 0.0;
	Real avg_density_err = 0.0;

	while ((!solved || (m_iterations < m_minIterations)) && m_iterations < m_maxIterations){
		solved = true;
		avg_density_err = 0.0;

		//////////////////////////////////////////////////////////////////////////
		// compute pressure accelerations
		for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
			FluidModel* fm = sim->getFluidModel(fluidModelIndex);
			const unsigned int numParticles = fm->numActiveParticles();

			for (int i = 0; i < numParticles; i++) {
				computePressureAccel(fluidModelIndex, i, fm->getDensity0(), m_simulationData.getPressureRho2Data());
			}
		}
		//////////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////////
		// update pressure values
		for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
			FluidModel* fm = sim->getFluidModel(fluidModelIndex);
			const unsigned int numParticles = fm->numActiveParticles();
			const Real density0 = fm->getDensity0();

			density_err = 0.0;
			for(int i = 0; i < numParticles; i++) {
				const Real sourceTerm = m_simulationData.getSourceTerm(fluidModelIndex, i);
				Real& pressure_i = m_simulationData.getPressure(fluidModelIndex, i);
				Real& diag_i = m_simulationData.getDiagElement(fluidModelIndex, i);
				const Real aij_pj = compute_aij_pj(fluidModelIndex, i);

				// update pressure
				if (abs(diag_i) > m_eps){
					pressure_i += (static_cast<Real>(0.5) / diag_i) * (sourceTerm - aij_pj);
					// pressure_i += (static_cast<Real>(0.5) / diag_i) * ((density0 - density_adv) * (static_cast<Real>(1.0) / h)  - aij_pj);
				}
				else {
					pressure_i = 0.0;
				}

				// clamp to zero
				pressure_i = max(pressure_i, 0.0f);
				density_err -= min((sourceTerm - aij_pj), 0.0f) * h;
				// density_err -= min((density0 - density_adv) * (static_cast<Real>(1.0) / h) - aij_pj, 0.0f) * h;
			}

			avg_density_err = density_err / numParticles;
			solved = solved && (avg_density_err <= (m_maxErrorV * 0.01 * density0));
		 }
		//////////////////////////////////////////////////////////////////////////
		///
		m_iterations++;
	}

	// final acceleration update
	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
		FluidModel* fm = sim->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = fm->numActiveParticles();

		for (int i = 0; i < numParticles; i++) {
			computePressureAccel(fluidModelIndex, i, fm->getDensity0(), m_simulationData.getPressureRho2Data());
		}
	}
}

void TimeStepDFSPHbubble::divergenceSolve()
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
		// compute pressure accelerations
		//////////////////////////////////////////////////////////////////////////
		for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
			FluidModel* fm = Simulation::getCurrent()->getFluidModel(fluidModelIndex);
			const unsigned int numParticles = fm->numActiveParticles();

			for (int i = 0; i < numParticles; i++) {
				computePressureAccel(fluidModelIndex, i, fm->getDensity0(), m_simulationData.getPressureRho2VData());
			}
		}

		//////////////////////////////////////////////////////////////////////////
		// update pressure values
		//////////////////////////////////////////////////////////////////////////
		for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
			FluidModel* fm = sim->getFluidModel(fluidModelIndex);
			const unsigned int numParticles = fm->numActiveParticles();
			Real divergence_error = 0.0;
			const Real density0 = fm->getDensity0();

			for (int i = 0; i < numParticles; i++){
				const Real sourceTermV = m_simulationData.getSourceTermDiv(fluidModelIndex, i);
				Real& pressureV_i = m_simulationData.getPressure_V(fluidModelIndex, i);
				const Real diag_i = m_simulationData.getDiagElement(fluidModelIndex, i);
				const Real aij_pj = compute_aij_pj(fluidModelIndex, i);

				// looking for neighbors in all point sets (incl. boundaries)
				unsigned int numNeighbors = 0;
				numNeighbors += sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i);

				// update pressure
				if(numNeighbors >= 7 && (abs(diag_i) > m_eps)){
					pressureV_i += (static_cast<Real>(0.5) / diag_i) * (sourceTermV - aij_pj);
				}
				else {
					pressureV_i = 0.0;
				}

				divergence_error -= min((sourceTermV - aij_pj), 0.0f) * h;
			}

			avg_divergence_error = divergence_error / numParticles;
			solved = solved && (avg_divergence_error <= (m_maxError * 0.01 * density0));
		}
		m_iterationsV++;
	}

	//////////////////////////////////////////////////////////////////////////
	// final pressure accel update
	//////////////////////////////////////////////////////////////////////////
	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
		FluidModel* fm = Simulation::getCurrent()->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = fm->numActiveParticles();

		for (int i = 0; i < numParticles; i++) {
			computePressureAccel(fluidModelIndex, i, fm->getDensity0(), m_simulationData.getPressureRho2VData());
		}
	}
}


void TimeStepDFSPHbubble::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
	m_iterations = 0;
	m_iterationsV = 0;
}

void TimeStepDFSPHbubble::performNeighborhoodSearch()
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

void TimeStepDFSPHbubble::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepDFSPHbubble::resize()
{
	m_simulationData.init();
}

/** Compute pressure accelerations using the current pressure values of the particles
 */
void TimeStepDFSPHbubble::computePressureAccel(const unsigned int fluidModelIndex, const unsigned int i, const Real density0, const std::vector<std::vector<Real>>& pressureList, const bool applyBoundaryForces)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int nFluids = sim->numberOfFluidModels();

	Vector3r& ap_i = m_simulationData.getPressureAccel(fluidModelIndex, i);
	ap_i.setZero();

	const Real& p_i = pressureList[fluidModelIndex][i];
	const Real density_i = model->getDensity(i);
	const Real dpi = p_i / (density_i * density_i);

	const Vector3r& xi = model->getPosition(i);

	
	forall_fluid_neighbors_in_same_phase (
		const Real& mass = model->getMass(neighborIndex);
		const Real& p_j = pressureList[fluidModelIndex][neighborIndex];
		const Real density_j = model->getDensity(neighborIndex);
		const Real dpj = p_j / (density_j * density_j);
		const Vector3r gradW = sim->gradW(xi - xj);

		ap_i -= mass * (dpi + dpj) * gradW;
	);

	//////////////////////////////////////////////////////////////////////////
	// Boundary Handling
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors(
			const Vector3r grad_p_j = /*psi*/ (bm_neighbor->getVolume(neighborIndex) * density0) * sim->gradW(xi - xj);
			const Vector3r a = (Real)1.0 * dpi * grad_p_j;
			ap_i -= a;
		);
	}
	//////////////////////////////////////////////////////////////////////////
}

void TimeStepDFSPHbubble::compute_aii(const unsigned int fluidModelIndex, const unsigned int i, const Real h){
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int nFluids = sim->numberOfFluidModels();

	const Real& rho_i = model->getDensity(i);
	const Real& density0_i = model->getDensity0();
	const Vector3r xi = model->getPosition(i);

	Real sumNormGradW = 0.0;
	Vector3r sumGrad = Vector3r::Zero();

	forall_fluid_neighbors_in_same_phase(
		Vector3r mgw = model->getMass(neighborIndex) * sim->gradW(xi - xj);
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
	m_simulationData.setDiagElement(fluidModelIndex, i, (-h / (rho_i * rho_i)) * sumNormGradW);
}

Real TimeStepDFSPHbubble::compute_aij_pj(const unsigned int fluidModelIndex, const unsigned int i)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const Real density0 = model->getDensity0();
	//////////////////////////////////////////////////////////////////////////
	// Compute A*p which is the change of the density when applying the
	// pressure forces.
	//////////////////////////////////////////////////////////////////////////
	const Vector3r& xi = model->getPosition(i);
	const Vector3r& ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
	Real aij_pj = 0.0;
	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors_in_same_phase(
		const Vector3r& aj = m_simulationData.getPressureAccel(fluidModelIndex, neighborIndex);
		aij_pj += model->getMass(neighborIndex) * (ai - aj).dot(sim->gradW(xi - xj));
		);
	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors(
			aij_pj += bm_neighbor->getVolume(neighborIndex) * density0 * ai.dot(sim->gradW(xi - xj));
		);
	}

	return aij_pj * h;
}


void TimeStepDFSPHbubble::computeConstantDensitySourceTerm (const unsigned int fluidModelIndex, const unsigned int i, const Real h){

	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	Real& sourceTerm = m_simulationData.getSourceTerm(fluidModelIndex, i);
	const Real density0 = model->getDensity0();
	const Real density = model->getDensity(i);
	const Vector3r& xi = model->getPosition(i);
	const Vector3r& vi = model->getVelocity(i);

	Real Drho_Dt = static_cast<Real>(0.0);

	forall_fluid_neighbors_in_same_phase (
        const Vector3r gradW = sim->gradW(xi - xj);
		Drho_Dt += model->getMass(i) * (vi - model->getVelocity(neighborIndex)).dot(gradW);
    );

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012){
		forall_boundary_neighbors(
			const Vector3r gradW = sim->gradW(xi - xj);
			Drho_Dt += (bm_neighbor->getVolume(neighborIndex) * density0) * vi.dot(gradW);
		);
	}
	//////////////////////////////////////////////////////////////////////////

	sourceTerm = 1.0 / h * (density0 - (density + h * Drho_Dt));
}

void TimeStepDFSPHbubble::computeDivergenceSourceTerm(const unsigned int fluidModelIndex, const unsigned int i, const Real h){
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int nFluids = sim->numberOfFluidModels();
	const Vector3r& xi = model->getPosition(i);

	Real& sourceTermDiv = m_simulationData.getSourceTermDiv(fluidModelIndex, i);
	Real denseDivVelocity = 0.0;

	forall_fluid_neighbors_in_same_phase(
		denseDivVelocity += model->getMass(i) * (model->getVelocity(i) - model->getVelocity(neighborIndex)).dot(sim->gradW(xi - xj));
	);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012) {
		forall_boundary_neighbors(
			denseDivVelocity += /* psi */ (bm_neighbor->getVolume(neighborIndex) * model->getDensity0()) * (model->getVelocity(i)).dot(sim->gradW(xi - xj));
		);
	}
	//////////////////////////////////////////////////////////////////////////

	sourceTermDiv = -denseDivVelocity;
}


// Two-way coupling
void TimeStepDFSPHbubble::computeDragForce(const unsigned int fluidModelIndex, const Real h){
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	unsigned int nFluids = sim->numberOfFluidModels();
	unsigned int numParticles = model->numActiveParticles();
	const Real h2 = h*h;

	for(int i = 0; i < numParticles; i++){
		const Real& density_i = model->getDensity(i);
		const Vector3r& vi = model->getVelocity(i);
		const Vector3r& xi = model->getPosition(i);

		Vector3r& acceleration = model->getAcceleration(i);

		Vector3r acc_drag = Vector3r::Zero();

		// drag constant of liquid is considered to be lower because influence of air onto liq might be lower than vice versa.
		Real dragConstant = m_dragConstantLiq;
		if (model->getId() != "Air"){
			dragConstant = m_dragConstantAir;
		}

		// Note: Does only work for Bubble framework if there is only one liquid and one gas!
		// So the Makro will only loop over the Air-phase for liquid particles and vice versa.
		forall_fluid_neighbors_in_different_phase(
			const Real m_j = fm_neighbor->getMass(neighborIndex);
			const Real density_j = fm_neighbor->getDensity(neighborIndex);
			const Vector3r& vj = fm_neighbor->getVelocity(neighborIndex);

			const Real pi_ij = max(0.0f, ((vi - vj).dot(xi - xj))/((xi - xj).norm() + m_eps * (h2)));
			// const Real pi_ij = max(0.0f, ((vi - vj).dot(xi - xj))/((xi - xj).squaredNorm() + m_eps * (h2)));

			acc_drag += m_j * ((dragConstant*h*m_speedSound)/(density_j + density_i)) * pi_ij * sim->gradW(xi - xj);
		);

		acceleration += acc_drag;
	}
}

void TimeStepDFSPHbubble::computeBouyancyForce(const unsigned int fluidModelIndex, const Real h){
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	unsigned int nFluids = sim->numberOfFluidModels();
	unsigned int numParticles = model->numActiveParticles();
	const Vector3r grav(sim->getVecValue<Real>(Simulation::GRAVITATION));

	for(int i = 0; i < numParticles; i++){
		Vector3r& acceleration = model->getAcceleration(i);

		Vector3r acc_bouyancy = Vector3r::Zero();

		// look for number of air particle neighbors
		int numNeighbors = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i);

		// equation (10)
		acc_bouyancy = m_minBouyancy * (m_kmax - (m_kmax - 1) * exp(-0.1 * numNeighbors)) * grav;
		acceleration -= acc_bouyancy;
	}
}

void TimeStepDFSPHbubble::computeCohesionForce(const unsigned int fluidModelIndex, const Real h) {
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	unsigned int nFluids = sim->numberOfFluidModels();
	unsigned int numParticles = model->numActiveParticles();

	for (int i = 0; i < numParticles; i++){
		Vector3r& acceleration = model->getAcceleration(i);
		const Vector3r& xi = model->getPosition(i);

		Vector3r acc_cohesion = Vector3r::Zero();

		forall_fluid_neighbors_in_same_phase(
			const Real densj = model->getDensity(neighborIndex);
			acc_cohesion += densj*(xi - xj);
		);

		acceleration -= m_cohesionConstant * acc_cohesion;
	}
}

void TimeStepDFSPHbubble::computeViscosityForce(const unsigned int fluidModelIndex, const unsigned i, const Real h) {
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = Simulation::getCurrent()->getFluidModel(fluidModelIndex);
	Vector3r& accel_i = model->getAcceleration(i);
	const Vector3r &xi = model->getPosition(i);

	Real factor = 0.01 / h;
	Vector3r sum = Vector3r::Zero();

	forall_fluid_neighbors_in_same_phase(
		sum += (model->getMass(neighborIndex) / model->getDensity(neighborIndex)) * (model->getVelocity(neighborIndex) - model->getVelocity(i)) * CubicKernel::W(xi - xj);
		);

	accel_i += factor * sum;
}

void TimeStepDFSPHbubble::computeSurfaceTensionForce(const unsigned int fluidModelIndex, const Real h){
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	unsigned int nFluids = sim->numberOfFluidModels();
	unsigned int numParticles = model->numActiveParticles();

	if(model->getId() == "Air") {
		LOG_ERR << "Surface tension force not implemented for air particles.";
		return;
	}

	for (int i = 0; i < numParticles; i++){
		Vector3r& acceleration = model->getAcceleration(i);
		const Real mi = model->getMass(i);
		const Vector3r xi = model->getPosition(i);

		Vector3r acc_surfaceTension = Vector3r::Zero();

		forall_fluid_neighbors_in_same_phase(
			const Real massj = model->getMass(neighborIndex);
			const Vector3r xij = xi-xj;

			acc_surfaceTension += massj * (xij) * sim->W(xij);
		);

		acceleration -= (m_surfaceTensionConstant * acc_surfaceTension) * (static_cast<Real>(1.0)/mi);
	}
}







