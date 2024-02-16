#include "TimeStepDFSPHbubbleOp.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include <iostream>
#include <random>
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"


using namespace SPH;
using namespace std;
using namespace GenParam;


int TimeStepDFSPHbubbleOp::SOLVER_ITERATIONS_V = -1;
int TimeStepDFSPHbubbleOp::MAX_ITERATIONS_V = -1;
int TimeStepDFSPHbubbleOp::MAX_ERROR_V = -1;
int TimeStepDFSPHbubbleOp::USE_DIVERGENCE_SOLVER = -1;
int TimeStepDFSPHbubbleOp::SOLVER_ITERATIONS_AIR = -1;
int TimeStepDFSPHbubbleOp::SOLVER_ITERATIONS_LIQ = -1;
int TimeStepDFSPHbubbleOp::SOLVER_ITERATIONS_V_AIR = -1;
int TimeStepDFSPHbubbleOp::SOLVER_ITERATIONS_V_LIQ = -1;
int TimeStepDFSPHbubbleOp::MAX_ERROR_AIR = -1;

int TimeStepDFSPHbubbleOp::USE_TRAPPED_AIR = -1;
int TimeStepDFSPHbubbleOp::USE_TRAPPED_AIR_OPTIMIZATION = -1;
int TimeStepDFSPHbubbleOp::VMIN_TRAPPED_AIR = -1;
int TimeStepDFSPHbubbleOp::VT_TRAPPED_AIR = -1;
int TimeStepDFSPHbubbleOp::VDIFF_THRESHOLD_MIN = -1;
int TimeStepDFSPHbubbleOp::VDIFF_THRESHOLD_MAX = -1;
int TimeStepDFSPHbubbleOp::MAX_AIR_PARTICLES_PER_STEP = -1;
int TimeStepDFSPHbubbleOp::EMIT_TIME_DISTANCE = -1;
int TimeStepDFSPHbubbleOp::DENSITY_RATIO_CAVITATION = -1;

int TimeStepDFSPHbubbleOp::TRAPPED_AIR_APPROACH = -1;
int TimeStepDFSPHbubbleOp::ENUM_TRAPPED_AIR_APPROACH_NONE = -1;
int TimeStepDFSPHbubbleOp::ENUM_TRAPPED_AIR_APPROACH_IHMSEN2011 = -1;
int TimeStepDFSPHbubbleOp::ENUM_TRAPPED_AIR_APPROACH_IHMSEN2012 = -1;
int TimeStepDFSPHbubbleOp::ENUM_TRAPPED_AIR_APPROACH_CAVITATION = -1;


TimeStepDFSPHbubbleOp::TimeStepDFSPHbubbleOp() :
	TimeStep(),
	m_simulationData()
{
	m_simulationData.init();
	m_counter = 0;
	m_iterationsV = 0;
	m_enableDivergenceSolver = true;
	m_enableTrappedAir = false;
	m_enableTrappedAirOptimization = false;
	m_maxIterationsV = 100;
	m_maxErrorV = static_cast<Real>(0.1);
	m_vMinTrappedAir = static_cast<Real>(9.0);
	m_vtTrappedAir = static_cast<Real>(0.3);
	m_nextEmitTime = static_cast<Real>(0.0);
	m_emitTimeDistance = static_cast<Real>(0.1);
	m_vDiffThresholdMin = static_cast<Real>(5.0);
	m_vDiffThresholdMax = static_cast<Real>(20.0);
	m_thresholdCavitationDensityRatio = static_cast<Real>(0.5);
	m_numberEmittedTrappedAirParticles = 0;
	m_iterationsAir = 0;
	m_iterationsLiq = 0;
	m_iterationsVair = 0;
	m_iterationsVliq = 0;
	m_trappedAirApproach = 1;
	m_maxErrorAir = 0.01;
	m_maxAirParticlesPerTimestep = 20;

	// add particle fields - then they can be used for the visualization and export
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->addField({ "factor", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getFactor(fluidModelIndex, i); } });
		model->addField({ "advected density", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getDensityAdv(fluidModelIndex, i); } });
		model->addField({ "p / rho^2", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureRho2(fluidModelIndex, i); }, true });
		model->addField({ "p_v / rho^2", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureRho2_V(fluidModelIndex, i); }, true });
		model->addField({ "pressure acceleration", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureAccel(fluidModelIndex, i)[0]; } });

		if (model->getId() == "Air")
			model->addField({ "lifetime", FieldType::Scalar, [this](const unsigned int i) -> Real* { return &m_simulationData.getLifetime(i); } });

	}
}

TimeStepDFSPHbubbleOp::~TimeStepDFSPHbubbleOp(void)
{
	// remove all particle fields
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		model->removeFieldByName("factor");
		model->removeFieldByName("advected density");
		model->removeFieldByName("p / rho^2");
		model->removeFieldByName("p_v / rho^2");
		model->removeFieldByName("pressure acceleration");

		if (model->getId() == "Air")
			model->removeFieldByName("lifetime");
	}

	LOG_INFO << "Sum of trapped Air particles: " << m_numberEmittedTrappedAirParticles;
}

void TimeStepDFSPHbubbleOp::initParameters()
{
	TimeStep::initParameters();

	SOLVER_ITERATIONS_AIR = createNumericParameter("iterationsAir", "Iterations density (air)", &m_iterationsAir);
	setGroup(SOLVER_ITERATIONS_AIR, "Simulation|DFSPH");
	setDescription(SOLVER_ITERATIONS_AIR, "Iterations required by the density solver for air model.");
	getParameter(SOLVER_ITERATIONS_AIR)->setReadOnly(true);

	SOLVER_ITERATIONS_LIQ = createNumericParameter("iterationsLiq", "Iterations density (liquid)", &m_iterationsLiq);
	setGroup(SOLVER_ITERATIONS_LIQ, "Simulation|DFSPH");
	setDescription(SOLVER_ITERATIONS_LIQ, "Iterations required by the density solver for liquid model.");
	getParameter(SOLVER_ITERATIONS_LIQ)->setReadOnly(true);

	SOLVER_ITERATIONS_V_AIR = createNumericParameter("iterationsVAir", "Iterations divergence (air)", &m_iterationsVair);
	setGroup(SOLVER_ITERATIONS_V_AIR, "Simulation|DFSPH");
	setDescription(SOLVER_ITERATIONS_V_AIR, "Iterations required by the divergence solver for air model.");
	getParameter(SOLVER_ITERATIONS_V_AIR)->setReadOnly(true);

	SOLVER_ITERATIONS_V_LIQ = createNumericParameter("iterationsVLiq", "Iterations divergence (liquid)", &m_iterationsVliq);
	setGroup(SOLVER_ITERATIONS_V_LIQ, "Simulation|DFSPH");
	setDescription(SOLVER_ITERATIONS_V_LIQ, "Iterations required by the divergence solver for liquid model.");
	getParameter(SOLVER_ITERATIONS_V_LIQ)->setReadOnly(true);

	SOLVER_ITERATIONS_V = createNumericParameter("iterationsV", "Iterations (divergence)", &m_iterationsV);
	setGroup(SOLVER_ITERATIONS_V, "Simulation|DFSPH");
	setDescription(SOLVER_ITERATIONS_V, "Iterations required by the divergence solver.");
	getParameter(SOLVER_ITERATIONS_V)->setReadOnly(true);

	MAX_ITERATIONS_V = createNumericParameter("maxIterationsV", "Max. iterations (divergence)", &m_maxIterationsV);
	setGroup(MAX_ITERATIONS_V, "Simulation|DFSPH");
	setDescription(MAX_ITERATIONS_V, "Maximal number of iterations of the divergence solver.");
	static_cast<NumericParameter<unsigned int>*>(getParameter(MAX_ITERATIONS_V))->setMinValue(1);

	MAX_ERROR_V = createNumericParameter("maxErrorV", "Max. divergence error(%)", &m_maxErrorV);
	setGroup(MAX_ERROR_V, "Simulation|DFSPH");
	setDescription(MAX_ERROR_V, "Maximal divergence error (%).");
	static_cast<RealParameter*>(getParameter(MAX_ERROR_V))->setMinValue(static_cast<Real>(1e-6));

	USE_DIVERGENCE_SOLVER = createBoolParameter("enableDivergenceSolver", "Enable divergence solver", &m_enableDivergenceSolver);
	setGroup(USE_DIVERGENCE_SOLVER, "Simulation|DFSPH");
	setDescription(USE_DIVERGENCE_SOLVER, "Turn divergence solver on/off.");

	TRAPPED_AIR_APPROACH = createEnumParameter("trappedAirApproach", "Trapped air approach", &m_trappedAirApproach);
	setGroup(TRAPPED_AIR_APPROACH, "Simulation|TrappedAir Extension");
	setDescription(TRAPPED_AIR_APPROACH, "Trapped air approach for air particle generation.");
	EnumParameter* enumParam = static_cast<EnumParameter*>(getParameter(TRAPPED_AIR_APPROACH));
	enumParam->addEnumValue("None", ENUM_TRAPPED_AIR_APPROACH_NONE);
	enumParam->addEnumValue("Ihmsen et al. 2011", ENUM_TRAPPED_AIR_APPROACH_IHMSEN2011);
	enumParam->addEnumValue("Ihmsen et al. 2012", ENUM_TRAPPED_AIR_APPROACH_IHMSEN2012);
	enumParam->addEnumValue("Cavitation Test", ENUM_TRAPPED_AIR_APPROACH_CAVITATION);

	USE_TRAPPED_AIR = createBoolParameter("enableTrappedAir", "Enable trapped air generation", &m_enableTrappedAir);
	setGroup(USE_TRAPPED_AIR, "Simulation|TrappedAir Extension");
	setDescription(USE_TRAPPED_AIR, "Turn trapped air generation on/off.");

	USE_TRAPPED_AIR_OPTIMIZATION = createBoolParameter("enableTrappedAirOptimization", "Enable trapped air optimization", &m_enableTrappedAirOptimization);
	setGroup(USE_TRAPPED_AIR_OPTIMIZATION, "Simulation|TrappedAir Extension");
	setDescription(USE_TRAPPED_AIR_OPTIMIZATION, "Turn trapped air optimization on/off - optimiation prevents emitation if distance to air particle is too small.");

	VMIN_TRAPPED_AIR = createNumericParameter("vMinTrappedAir", "Min. velocity trapped air", &m_vMinTrappedAir);
	setGroup(VMIN_TRAPPED_AIR, "Simulation|TrappedAir Extension");
	setDescription(VMIN_TRAPPED_AIR, "Minimal velocity of liquid particle to generate an air particle.");

	VT_TRAPPED_AIR = createNumericParameter("vtTrappedAir", "Ihmsen 2011: Threshold velocity difference", &m_vtTrappedAir);
	setGroup(VT_TRAPPED_AIR, "Simulation|TrappedAir Extension");
	setDescription(VT_TRAPPED_AIR, "Threshold velocity difference between liquid particle and its liquid neighbors for air particle generation.");
	
	VDIFF_THRESHOLD_MIN = createNumericParameter("vDiffThresholdMin", "Ihmsen 2012: Min. threshold velocity difference", &m_vDiffThresholdMin);
	setGroup(VDIFF_THRESHOLD_MIN, "Simulation|TrappedAir Extension");
	setDescription(VDIFF_THRESHOLD_MIN, "Min. threshold velocity difference between liquid particle and its liquid neighbors for air particle generation.");

	VDIFF_THRESHOLD_MAX = createNumericParameter("vDiffThresholdMax", "Ihmsen 2012: Max. threshold velocity difference", &m_vDiffThresholdMax);
	setGroup(VDIFF_THRESHOLD_MAX, "Simulation|TrappedAir Extension");
	setDescription(VDIFF_THRESHOLD_MAX, "Max. threshold velocity difference between liquid particle and its liquid neighbors for air particle generation.");

	MAX_ERROR_AIR = createNumericParameter("maxErrorAir", "Max. density error(%) for air phase", &m_maxErrorAir);
	setGroup(MAX_ERROR_AIR, "Simulation|DFSPHbubble");
	setDescription(MAX_ERROR_AIR, "Maximal density error (%) for air phase.");

	MAX_AIR_PARTICLES_PER_STEP = createNumericParameter("maxAirParticlesPerStep", "Max. num air particles generated per step", &m_maxAirParticlesPerTimestep);
	setGroup(MAX_AIR_PARTICLES_PER_STEP, "Simulation|TrappedAir Extension");
	setDescription(MAX_AIR_PARTICLES_PER_STEP, "Max. number of air particles that can be generated per emitting time step.");

	EMIT_TIME_DISTANCE = createNumericParameter("emitTimeDistance", "Temporal distance between trapped Air generation steps (seconds)", &m_emitTimeDistance);
	setGroup(EMIT_TIME_DISTANCE, "Simulation|TrappedAir Extension");
	setDescription(EMIT_TIME_DISTANCE, "Temporal distance between trapped Air generation steps (seconds)");

	DENSITY_RATIO_CAVITATION = createNumericParameter("densityRatioCavitation", "Threshold density ratio for cavitation trapped air approach", &m_thresholdCavitationDensityRatio);
	setGroup(DENSITY_RATIO_CAVITATION, "Simulation|TrappedAir Extension");
	setDescription(DENSITY_RATIO_CAVITATION, "Threshold density ratio for cavitation trapped air approach.");
}

void TimeStepDFSPHbubbleOp::step()
{
	Simulation *sim = Simulation::getCurrent();
	TimeManager *tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();
	const unsigned int nModels = sim->numberOfFluidModels();

	//////////////////////////////////////////////////////////////////////////
	// search the neighbors for all particles
	//////////////////////////////////////////////////////////////////////////
	performNeighborhoodSearch();

#ifdef USE_PERFORMANCE_OPTIMIZATION
	//////////////////////////////////////////////////////////////////////////
	// pre-compute the values V_j * grad W_ij for all neighbors
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("precomputeValues")
	precomputeValues();
	STOP_TIMING_AVG
#endif

	//////////////////////////////////////////////////////////////////////////
	// compute volume/density maps boundary contribution
	//////////////////////////////////////////////////////////////////////////
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		computeVolumeAndBoundaryX();
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		computeDensityAndGradient();

	//////////////////////////////////////////////////////////////////////////
	// compute densities
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
		computeDensitiesForSamePhaseAVX(fluidModelIndex);

	//////////////////////////////////////////////////////////////////////////
	// Compute the factor alpha_i for all particles i
	// using the equation (11) in [BK17]
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("computeDFSPHFactor");
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++){
		computeDFSPHFactor(fluidModelIndex);
	}
	STOP_TIMING_AVG;

	//////////////////////////////////////////////////////////////////////////
	// Perform divergence solve (see Algorithm 2 in [BK17])
	//////////////////////////////////////////////////////////////////////////
	if (m_enableDivergenceSolver)
	{
		START_TIMING("divergenceSolve");
		divergenceSolve();
		STOP_TIMING_AVG
	}
	else {
		m_iterationsV = 0;
		m_iterationsVair = 0;
	}

	//////////////////////////////////////////////////////////////////////////
	// Reset accelerations and add gravity
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++) {
		clearAccelerations(fluidModelIndex);
	}


	//////////////////////////////////////////////////////////////////////////
	// Track state of fluid models
	//////////////////////////////////////////////////////////////////////////
	int liquidModelIndex = -1;
	int airModelIndex = -1;

	if(nModels == 2){
		liquidModelIndex = sim->getFluidModel(0)->getId() == "Liquid" ? sim->getFluidModel(0)->getPointSetIndex() : sim->getFluidModel(1)->getPointSetIndex();
		airModelIndex = sim->getFluidModel(0)->getId() == "Air" ? sim->getFluidModel(0)->getPointSetIndex() : sim->getFluidModel(1)->getPointSetIndex();
	}
	else if(nModels == 1){
		if(sim->getFluidModel(0)->getId() == "Liquid"){
			liquidModelIndex = 0;
		}
		else{
			airModelIndex = 0;
		}
	}
	else{
		LOG_ERR << "No or too many fluid-models not supported yet. Number fluid-models: " << nModels;
		return;
	}

	unsigned int nLiquidParticles = 0;
	if(liquidModelIndex > -1){
		nLiquidParticles = sim->getFluidModel(liquidModelIndex)->numActiveParticles();
	}

	unsigned int nAirParticles = 0;
	if(airModelIndex > -1){
		nAirParticles = sim->getFluidModel(airModelIndex)->numActiveParticles();
	}

	//////////////////////////////////////////////////////////////////////////
	// Non-pressure forces
	//////////////////////////////////////////////////////////////////////////
	sim->computeNonPressureForces();

	// liquid pase exclusive non-pressure forces
	// if(nLiquidParticles > 0){
	// 	//////////////////////////////////////////////////////////////////////////
	// 	// Viscosity XSPH -> only liquid model
	// 	//////////////////////////////////////////////////////////////////////////
	// 	#pragma omp parallel default(shared)
	// 	{
	// 		#pragma omp for schedule(static)
	// 		for(int i = 0; i < nLiquidParticles; i++){
	// 			computeViscosityForce(liquidModelIndex, i, h);
	// 		}
	// 	}
	// 
	// 	//////////////////////////////////////////////////////////////////////////
	// 	// Surface Tension Force -> for liquid model
	// 	//////////////////////////////////////////////////////////////////////////
	// 	computeSurfaceTensionForce(liquidModelIndex, h);
	// }

	//////////////////////////////////////////////////////////////////////////
	// Non-pressure Forces introduced by Bubble-paper Ihmsen et al.
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("bubbleForces");
	if(nLiquidParticles > 0){
		sim->getFluidModel(liquidModelIndex)->computeBubbleForces();
	}
	if(nAirParticles > 0){
		//////////////////////////////////////////////////////////////////////////
		// Compute onSurface state and lifetime updates for air particles
		//////////////////////////////////////////////////////////////////////////
		computeOnSurfaceAir();

		sim->getFluidModel(airModelIndex)->computeBubbleForces();
	}
	STOP_TIMING_AVG;

	//////////////////////////////////////////////////////////////////////////
	// Update the time step size, e.g. by using a CFL condition
	//////////////////////////////////////////////////////////////////////////
	sim->updateTimeStepSize();

	//////////////////////////////////////////////////////////////////////////
	// compute new velocities only considering non-pressure forces
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int m = 0; m < nModels; m++)
	{
		FluidModel *fm = sim->getFluidModel(m);
		const unsigned int numParticles = fm->numActiveParticles();
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				if (fm->getParticleState(i) == ParticleState::Active)
				{
					Vector3r &vel = fm->getVelocity(i);
					vel += h * fm->getAcceleration(i);
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Perform constant density solve (see Algorithm 3 in [BK17])
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("pressureSolve");
	pressureSolve();
	STOP_TIMING_AVG;

	//////////////////////////////////////////////////////////////////////////
	// compute final positions
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int m = 0; m < nModels; m++)
	{
		FluidModel *fm = sim->getFluidModel(m);
		unsigned int deletedParticles = 0;
		const unsigned int numParticles = fm->numActiveParticles();
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)
			for (int i = 0; i < ((int)numParticles - deletedParticles); i++)
			{
				if (fm->getParticleState(i) == ParticleState::Active)
				{
					Vector3r &xi = fm->getPosition(i);
					const Vector3r &vi = fm->getVelocity(i);
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
	// air particle generation: Trapped Air
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("trappedAirGeneration");
	if (m_enableTrappedAir && nLiquidParticles > 0 && liquidModelIndex > -1) {
		FluidModel* liquid = sim->getFluidModel(liquidModelIndex);
		const Real t = tm->getTime();

		if (t >= m_nextEmitTime) {
			// loop over liquid particle
			unsigned int emittedParticles = 0;
			std::vector<unsigned int> indicesGen;
			indicesGen.clear();

			if (m_trappedAirApproach == ENUM_TRAPPED_AIR_APPROACH_IHMSEN2011) {
				for (unsigned int i = 0; i < nLiquidParticles; i++) {
					if (emittedParticles < m_maxAirParticlesPerTimestep && liquid->getParticleState(i) == ParticleState::Active) {
						trappedAirIhmsen2011(liquidModelIndex, i, emittedParticles, indicesGen);
					}
				}
			}
			else if (m_trappedAirApproach == ENUM_TRAPPED_AIR_APPROACH_IHMSEN2012) {
				trappedAirIhmsen2012(emittedParticles, indicesGen);
			}
			else if (m_trappedAirApproach == ENUM_TRAPPED_AIR_APPROACH_CAVITATION) {
				for (unsigned int i = 0; i < nLiquidParticles; i++) {
					if (emittedParticles < m_maxAirParticlesPerTimestep && liquid->getParticleState(i) == ParticleState::Active) {
						trappedAirCavitation(liquidModelIndex, i, emittedParticles, indicesGen);
					}
				}
			}

			if (emittedParticles != 0)
			{
				FluidModel* airModel = sim->getFluidModel(airModelIndex);
				airModel->setNumActiveParticles(nAirParticles + emittedParticles);
				//sim->emittedParticles(airModel, airModel->numActiveParticles() - emittedParticles);
				this->emittedParticles(airModel, airModel->numActiveParticles() - emittedParticles);
				sim->getNeighborhoodSearch()->resize_point_set(airModel->getPointSetIndex(), &airModel->getPosition(0)[0], airModel->numActiveParticles());

				m_nextEmitTime += m_emitTimeDistance;
			}

			indicesGen.clear();
		}

	}
	STOP_TIMING_AVG;

	//////////////////////////////////////////////////////////////////////////
	// Compute new time
	//////////////////////////////////////////////////////////////////////////
	tm->setTime (tm->getTime() + h);
}




//////////////////////////////////////////////////////////////////////////
// emit air particles based on the velocity field of the liquid -> Ihmsen et al. 2011
//////////////////////////////////////////////////////////////////////////
void TimeStepDFSPHbubbleOp::trappedAirIhmsen2011(const unsigned int liquidModelIndex, const unsigned int i, unsigned int& emittedParticles, std::vector<unsigned int>& indicesGen){
	Simulation* sim = Simulation::getCurrent();
	FluidModel* airModel = sim->getFluidModel(1)->getId() == "Air" ? sim->getFluidModel(1) : sim->getFluidModel(0);
	FluidModel* model = sim->getFluidModel(liquidModelIndex); 	// liquid model -> naming convention for makros
	const unsigned int numLiqParticles = model->numActiveParticles();
	const unsigned int fluidModelIndex = liquidModelIndex;

	if((airModel->numActiveParticles() >= airModel->numParticles()))
		return;

	const Vector3r xi = model->getPosition(i);
	const Vector3r vi = model->getVelocity(i);

	if(vi.squaredNorm() < m_vMinTrappedAir){
		return;
	}

	//////////////////////////////////////////////////////////////////////////
	// extension AK: check if there is any air particle too close to the current liquid particle
	if (m_enableTrappedAirOptimization) {
		const Real radius = sim->getParticleRadius();
		const Real diam = 2 * radius;
		const unsigned int nFluids = sim->numberOfFluidModels();

		volatile bool isTooClose = false;
		// loop over all air particles in this specific case
		forall_fluid_neighbors_in_different_phase(
			if (isTooClose) {
				continue;
			}

		const Vector3r xij = xi - xj;
		if (xij.squaredNorm() < (diam * diam)) {
			isTooClose = true;
		}
		);

		// check newly generated air particle positions
		for (unsigned int j : indicesGen) {
			if (isTooClose) {
				continue;
			}

			const Vector3r xij = xi - model->getPosition(j);
			if (xij.squaredNorm() < (diam * diam)) {
				isTooClose = true;
			}
		}

		if (isTooClose) {
			// neighbor air particle too close -> skip this liquid particle
			return;
		}
	}
	//////////////////////////////////////////////////////////////////////////

	// compute v_diff
	Vector3r v_diff = Vector3r::Zero();
	
	// loop over liquid neighbors
	forall_fluid_neighbors_in_same_phase(
		const Vector3r& vj = model->getVelocity(neighborIndex);
		const Real& density_j = model->getDensity(neighborIndex);
		const Vector3r v_ij = vi - vj;
	
		v_diff += (1.0 / density_j) * v_ij * sim->W(xi - xj);
	);
	
	v_diff *= model->getMass(i);
	Real vDiffNorm = v_diff.norm();
	
	if (vDiffNorm < m_vtTrappedAir){
		return;
	}
	
	// get number of air neighbors of liquid particle i
	unsigned int numAirNeighbors = 0;
	numAirNeighbors += sim->numberOfNeighbors(fluidModelIndex, airModel->getPointSetIndex(), i);
	
	if (numAirNeighbors >= (vDiffNorm / m_vtTrappedAir))
		return;

	// emit air particle
	emitAirParticleFromVelocityField(emittedParticles, vi, xi);
	indicesGen.push_back(i);
}


// trapped air generation Ihmsen et al. 2012: "Unified spray, foam and air bubbles for particle-based fluids"
void TimeStepDFSPHbubbleOp::trappedAirIhmsen2012(unsigned int& emittedParticles, std::vector<unsigned int>& indicesGen)
{
	const Real v_min = m_vMinTrappedAir; // min. velocity for liquid particle to generate an air particle, original: 9.0
	const Real tMin = m_vDiffThresholdMin; // original 5.0
	const Real tMax = m_vDiffThresholdMax; // original 20.0

	auto clamp = [&tMin, &tMax](const Real& v_diff){
		return (min(v_diff, tMax) - min(v_diff, tMin)) / (tMax - tMin);
	};


	Simulation* sim = Simulation::getCurrent();
	const Real& supportRadius = sim->getSupportRadius();
	FluidModel* airModel = sim->getFluidModel(1)->getId() == "Air" ? sim->getFluidModel(1) : sim->getFluidModel(0);
	FluidModel* model = sim->getFluidModel(1)->getId() == "Air" ? sim->getFluidModel(0) : sim->getFluidModel(1); // liquid model -> naming convention

	const unsigned int fluidModelIndex = model->getPointSetIndex();
	const unsigned int numLiqParticles = model->numActiveParticles();

	auto customW = [&supportRadius](const Vector3r& dist) -> Real{
		const Real distNorm = dist.norm();
		if(distNorm <= supportRadius){
			return (1.0 - ((1.0/supportRadius)*distNorm));
		}
		return 0.0;
	};


	// reached limit of air particle in simulation
	if ((airModel->numActiveParticles() >= airModel->numParticles()))
		return;

	// loop over liquid particles
	for(int i = 0; i < numLiqParticles; i++){

		if (emittedParticles >= m_maxAirParticlesPerTimestep) {
			return;
		}

		if(model->getParticleState(i) != ParticleState::Active){
			continue;
        }

		const Vector3r& xi = model->getPosition(i);
		const Vector3r& vi = model->getVelocity(i);

		// minimal required velocity for air emitting
		if(vi.squaredNorm() < v_min)
			continue;

		//////////////////////////////////////////////////////////////////////////
		// extension AK: check if there is any air particle too close to the current liquid particle
		if (m_enableTrappedAirOptimization) {
			const Real radius = sim->getParticleRadius();
			const Real diam = 2 * radius;
			const unsigned int nFluids = sim->numberOfFluidModels();

			volatile bool isTooClose = false;
			// loop over all air particles in this specific case
			forall_fluid_neighbors_in_different_phase(
				if (isTooClose) {
					continue;
				}

				const Vector3r xij = xi - xj;
				if (xij.squaredNorm() < (diam * diam)) {
					isTooClose = true;
				}
			);

			for (unsigned int k : indicesGen) {
				if (isTooClose) {
                    continue;
                }

                const Vector3r xik = xi - model->getPosition(k);
				if (xik.squaredNorm() < (diam * diam)) {
                    isTooClose = true;
                }
            }

			if (isTooClose) {
				// neighbor air particle too close -> skip this liquid particle
				continue;
			}
		}
		//////////////////////////////////////////////////////////////////////////

		Real probability = 0;
		Real v_diff = 0;

		forall_fluid_neighbors_in_same_phase(
			const Vector3r & vj = model->getVelocity(neighborIndex);

			const Vector3r vij = vi - vj;
			const Vector3r xij = xi - xj;
			const Vector3r vij_normalized = vij.normalized();
			const Vector3r xij_normalized = xij.normalized();

			v_diff += vij.norm() * (1 - vij_normalized.dot(xij_normalized)) * customW(xij);
		);

		probability = clamp(v_diff);

		// generate a particle if probability is greater than threshold
		if(probability < 0.5)
			continue;

		// emit air particle
		emitAirParticleFromVelocityField(emittedParticles, vi, xi);
		indicesGen.push_back(i);
	}
}

// trapped air via cavitation phenomenon in low density areas. See Lugli and Zerbetto 2013: "An introduction to bubble dynamics"
void SPH::TimeStepDFSPHbubbleOp::trappedAirCavitation(const unsigned int liquidModelIndex, const unsigned int i, unsigned int& emittedParticles, std::vector<unsigned int>& indicesGen)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* airModel = sim->getFluidModel(1)->getId() == "Air" ? sim->getFluidModel(1) : sim->getFluidModel(0);
	FluidModel* model = sim->getFluidModel(liquidModelIndex); 	// liquid model -> naming convention for makros
	const unsigned int numLiqParticles = model->numActiveParticles();
	const unsigned int fluidModelIndex = liquidModelIndex;
	const Vector3r grav(sim->getVecValue<Real>(Simulation::GRAVITATION));


	// requirements
	const Real density0 = model->getDensity0();
	const Real density_i = model->getDensity(i);

	const Real densRatio = density_i / density0;
	
	if (densRatio > m_thresholdCavitationDensityRatio) {
		return; // cavitation in areas with density significantly lower than rest density
	}

	const Vector3r xi = model->getPosition(i);

	const unsigned int th_onSurface = static_cast<unsigned int>(3); // required nbr of liquid particles of curr particle to be considered "onSurface"
	unsigned int onSurfaceCount = 0; // keeping track how many liquid particles are above curr. particle
	bool onSurface = true;

	// additionally check if particle on surface (see Ihmsen et al. 2011)
	forall_fluid_neighbors_in_same_phase(
		if (!onSurface) {
			continue;
		}
		if ((xj - xi).dot(grav) < 0) {
			onSurfaceCount++;
		}
		if (onSurfaceCount >= th_onSurface) {
			onSurface = false;
		}
	);

	if (onSurface) {
		return; // no cavitation on surface of a fluid
	}

	const Vector3r vi = model->getVelocity(i);

	emitAirParticleFromVelocityField(emittedParticles, vi, xi);
	indicesGen.push_back(i);
}

// compute state of the air particles: on surface or inside liquid
void TimeStepDFSPHbubbleOp::computeOnSurfaceAir(){
	Simulation *sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(0)->getId() == "Air" ? sim->getFluidModel(0) : sim->getFluidModel(1); // air model
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	const Vector3r grav(sim->getVecValue<Real>(Simulation::GRAVITATION));
	TimeManager* tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	// compute onSurface flag and updates of lifetime (multithreaded)
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for(int i = 0; i < numParticles; i++)
		{
			if(model->getParticleState(i) == ParticleState::Disabled){
				continue;
			}

			const Real density_i = model->getDensity(i);
			const Vector3r& xi = model->getPosition(i);

			// This condition seems to be error prone. If an air particle gets "trapped" it might not have any air-neighbors and would be falsly identified as "on the surface"
			// if(density_i > m_onSurfaceThresholdDensity){
			// 	m_onSurface[i] = 1;
			// }

			// look for number of liquid particle neighbors of the air particle
			// int numLiqNeighbors = sim->numberOfNeighbors(m_model->getPointSetIndex(), fluidModelIndex, i);

			volatile bool onSurface = true;
			int liquidNeighbors = 0;
			// looping over liquid neighbors
			forall_fluid_neighbors_in_different_phase(
				if(!onSurface){
					continue;
				}
				if((xj-xi).dot(grav) < 0){
					onSurface = false;
				}
				liquidNeighbors++;
			);

			if (onSurface && m_simulationData.getOnSurface(i) == 0) {
				Vector3r& vel = model->getVelocity(i);
				vel *= 0.1;
            }
			sim->getSupportRadius();
			m_simulationData.getOnSurface(i) = onSurface;
			Real& lifetime_i = m_simulationData.getLifetime(i);

			// lifetime update
			if(onSurface){
				if (liquidNeighbors <= 4 && lifetime_i > h) {
					lifetime_i = 2*h;
				}
				
				lifetime_i -= h;
			}
		}
	}

	// foam deletion check
	for(int i = 0; i < numParticles; i++){
		if(model->getParticleState(i) == ParticleState::Disabled){
			continue;
		}

		Real& lifetime_i = m_simulationData.getLifetime(i);

		// all particles of a bubble should disperse at once.
		forall_fluid_neighbors_in_same_phase(
			lifetime_i = std::min(lifetime_i, m_simulationData.getLifetime(neighborIndex));
		);

		// Disable an air particle at the end of its lifetime
		if(lifetime_i <= 0.0){
			model->getVelocity(i) = Vector3r::Zero(); // prevent that velocity influences the CFL condition
			model->getPosition(i) = Vector3r(1000 + i, 1000, 1000);

			model->setParticleState(i, ParticleState::Disabled);
			// -> clean-up in TimeStep
			// TODO: Clean-up can now take place here
			m_simulationData.setLifetime(i, 3.0);
			m_simulationData.getOnSurface(i) = 0;
		}
	}

}

void TimeStepDFSPHbubbleOp::pressureSolve()
{
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real h2 = h*h;
	const Real invH = static_cast<Real>(1.0) / h;
	const Real invH2 = static_cast<Real>(1.0) / h2;
	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();

	//////////////////////////////////////////////////////////////////////////
	// Compute rho_adv
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const Real density0 = model->getDensity0();
		const int numParticles = (int)model->numActiveParticles();
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
				if(model->getParticleState(i) == ParticleState::Disabled){
					continue;
				}
				//////////////////////////////////////////////////////////////////////////
				// Compute rho_adv,i^(0) (see equation in Section 3.3 in [BK17])
				// using the velocities after the non-pressure forces were applied.
				//////////////////////////////////////////////////////////////////////////
				computeDensityAdv(fluidModelIndex, i, h, density0);

				//////////////////////////////////////////////////////////////////////////
				// In the end of Section 3.3 [BK17] we have to multiply the density 
				// error with the factor alpha_i divided by h^2. Hence, we multiply 
				// the factor directly by 1/h^2 here.
				//////////////////////////////////////////////////////////////////////////
				m_simulationData.getFactor(fluidModelIndex, i) *= invH2;

				//////////////////////////////////////////////////////////////////////////
				// For the warm start we use 0.5 times the old pressure value.
				// Note: We divide the value by h^2 since we multiplied it by h^2 at the end of 
				// the last time step to make it independent of the time step size.
				//////////////////////////////////////////////////////////////////////////
#ifdef USE_WARMSTART
				if (m_simulationData.getDensityAdv(fluidModelIndex, i) > 1.0)
					m_simulationData.getPressureRho2(fluidModelIndex, i) = static_cast<Real>(0.5) * min(m_simulationData.getPressureRho2(fluidModelIndex, i), static_cast<Real>(0.00025)) * invH2;
				else 
					m_simulationData.getPressureRho2(fluidModelIndex, i) = 0.0;
#else 
				//////////////////////////////////////////////////////////////////////////
				// If we don't use a warm start, we directly compute a pressure value
				// by multiplying the density error with the factor.
				//////////////////////////////////////////////////////////////////////////
				//m_simulationData.getPressureRho2(fluidModelIndex, i) = 0.0;
				const Real s_i = static_cast<Real>(1.0) - m_simulationData.getDensityAdv(fluidModelIndex, i);
				const Real residuum = min(s_i, static_cast<Real>(0.0));     // r = b - A*p
				m_simulationData.getPressureRho2(fluidModelIndex, i) = -residuum * m_simulationData.getFactor(fluidModelIndex, i);
#endif
			}
		}
	}


	//////////////////////////////////////////////////////////////////////////
	// Start solver
	//////////////////////////////////////////////////////////////////////////
	Real avg_density_err = 0.0;
	bool chk = false;
	unsigned int iterations[2] = { 0,0 };
	assert(nFluids == 2);
	m_iterations = 0;

	for (unsigned int i = 0; i < nFluids; i++) {
		FluidModel* model = sim->getFluidModel(i);
		const Real density0 = model->getDensity0();
		iterations[i] = 0;
		Real maxErr = model->getId() == "Air" ? m_maxErrorAir : m_maxError;

		//////////////////////////////////////////////////////////////////////////
		// Perform solver iterations
		//////////////////////////////////////////////////////////////////////////
		while ((!chk || (iterations[i] < m_minIterations)) && (iterations[i] < m_maxIterations)) {
			chk = true;
			avg_density_err = 0.0;

			pressureSolveIteration(i, avg_density_err);

			// Maximal allowed density fluctuation
			const Real eta = maxErr * static_cast<Real>(0.01) * density0;  // maxError is given in percent
			chk = avg_density_err <= eta;

			iterations[i]++;
		}
	}

	m_iterations = iterations[0] + iterations[1];
	m_iterationsAir = sim->getFluidModel(0)->getId() == "Air" ? iterations[0] : iterations[1];
	m_iterationsLiq = sim->getFluidModel(0)->getId() == "Liquid" ? iterations[0] : iterations[1];

	INCREASE_COUNTER("DFSPH - iterations", static_cast<Real>(m_iterations));

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const int numParticles = (int)model->numActiveParticles();
		const Real density0 = model->getDensity0();
		
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
				//////////////////////////////////////////////////////////////////////////
				// Time integration of the pressure accelerations to get new velocities
				//////////////////////////////////////////////////////////////////////////
				computePressureAccel(fluidModelIndex, i, density0, m_simulationData.getPressureRho2Data(), true);
				model->getVelocity(i) += h * m_simulationData.getPressureAccel(fluidModelIndex, i);
			}
		}
	}
#ifdef USE_WARMSTART
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		const int numParticles = (int)model->numActiveParticles();
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
				//////////////////////////////////////////////////////////////////////////
				// Multiply by h^2, the time step size has to be removed 
				// to make the pressure value independent 
				// of the time step size
				//////////////////////////////////////////////////////////////////////////
				m_simulationData.getPressureRho2(fluidModelIndex, i) *= h2;
			}		
		}
	}
#endif
}

void TimeStepDFSPHbubbleOp::divergenceSolve()
{
	//////////////////////////////////////////////////////////////////////////
	// Init parameters
	//////////////////////////////////////////////////////////////////////////

	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real invH = static_cast<Real>(1.0) / h;
	Simulation *sim = Simulation::getCurrent();
	const unsigned int maxIter = m_maxIterationsV;
	const Real maxError = m_maxErrorV;
	const unsigned int nFluids = sim->numberOfFluidModels();

	//////////////////////////////////////////////////////////////////////////
	// Compute divergence of velocity field
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const int numParticles = (int)model->numActiveParticles();

		#pragma omp parallel default(shared)
		{		
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
				if(model->getParticleState(i) == ParticleState::Disabled){
					continue;
				}
				//////////////////////////////////////////////////////////////////////////
				// Compute rho_adv,i^(0) (see equation (9) in Section 3.2 [BK17])
				// using the velocities after the non-pressure forces were applied.
				//////////////////////////////////////////////////////////////////////////
				computeDensityChange(fluidModelIndex, i, h);

				Real densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
				densityAdv = max(densityAdv, static_cast<Real>(0.0));

				unsigned int numNeighbors = 0;
				// only consider neighbors of the own phase
				numNeighbors += sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i);

				// in case of particle deficiency do not perform a divergence solve
				if (!sim->is2DSimulation())
				{
					if (numNeighbors < 20)
						densityAdv = 0.0;
				}
				else
				{
					if (numNeighbors < 7)
						densityAdv = 0.0;
				}
				
				//////////////////////////////////////////////////////////////////////////
				// In equation (11) [BK17] we have to multiply the divergence 
				// error with the factor divided by h. Hence, we multiply the factor
				// directly by 1/h here.
				//////////////////////////////////////////////////////////////////////////
				m_simulationData.getFactor(fluidModelIndex, i) *= invH;

				//////////////////////////////////////////////////////////////////////////
				// For the warm start we use 0.5 times the old pressure value.
				// Divide the value by h. We multiplied it by h at the end of 
				// the last time step to make it independent of the time step size.
				//////////////////////////////////////////////////////////////////////////
#ifdef USE_WARMSTART_V
				if (densityAdv > 0.0)
					m_simulationData.getPressureRho2_V(fluidModelIndex, i) = static_cast<Real>(0.5) * min(m_simulationData.getPressureRho2_V(fluidModelIndex, i), static_cast<Real>(0.5)) * invH;
				else
					m_simulationData.getPressureRho2_V(fluidModelIndex, i) = 0.0;
#else 
				//////////////////////////////////////////////////////////////////////////
				// If we don't use a warm start, directly compute a pressure value
				// by multiplying the divergence error with the factor.
				//////////////////////////////////////////////////////////////////////////
				m_simulationData.getPressureRho2_V(fluidModelIndex, i) = densityAdv * m_simulationData.getFactor(fluidModelIndex, i);
#endif
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Start solver
	//////////////////////////////////////////////////////////////////////////	
	Real avg_density_err = 0.0;
	bool chk = false;

	unsigned int iterationsV[2] = { 0,0 };
	assert(nFluids == 2);
	m_iterationsV = 0;

	for (unsigned int i = 0; i < nFluids; i++) {
		iterationsV[i] = 0;
		FluidModel* model = sim->getFluidModel(i);
		const Real density0 = model->getDensity0();

		//////////////////////////////////////////////////////////////////////////
		// Perform solver iterations
		//////////////////////////////////////////////////////////////////////////
		while ((!chk || (iterationsV[i] < 1)) && (iterationsV[i] < maxIter))
		{
			chk = true;
			avg_density_err = 0.0;

			divergenceSolveIteration(i, avg_density_err);

			// Maximal allowed density fluctuation
			// use maximal density error divided by time step size
			const Real eta = (static_cast<Real>(1.0) / h) * maxError * static_cast<Real>(0.01) * density0;  // maxError is given in percent
			chk = avg_density_err <= eta;

			iterationsV[i]++;
		}
	}

	m_iterationsV = iterationsV[0] + iterationsV[1];
	m_iterationsVair = sim->getFluidModel(0)->getId() == "Air" ? iterationsV[0] : iterationsV[1];
	m_iterationsVliq = sim->getFluidModel(1)->getId() == "Liquid" ? iterationsV[1] : iterationsV[0];


	INCREASE_COUNTER("DFSPH - iterationsV", static_cast<Real>(m_iterationsV));

	
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel *model = sim->getFluidModel(fluidModelIndex);
		const int numParticles = (int)model->numActiveParticles();
		const Real density0 = model->getDensity0();

		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
				//////////////////////////////////////////////////////////////////////////
				// Time integration of the pressure accelerations
				//////////////////////////////////////////////////////////////////////////
				computePressureAccel(fluidModelIndex, i, density0, m_simulationData.getPressureRho2VData(), true);
				model->getVelocity(i) += h * m_simulationData.getPressureAccel(fluidModelIndex, i);

				m_simulationData.getFactor(fluidModelIndex, i) *= h;
			}
		}
	}
#ifdef USE_WARMSTART_V
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		const int numParticles = (int)model->numActiveParticles();
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < numParticles; i++)
			{
				//////////////////////////////////////////////////////////////////////////
				// Multiply by h, the time step size has to be removed 
				// to make the pressure value independent 
				// of the time step size
				//////////////////////////////////////////////////////////////////////////		
				m_simulationData.getPressureRho2_V(fluidModelIndex, i) *= h;
			}
		}
	}
#endif
}


void TimeStepDFSPHbubbleOp::pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const Real density0 = model->getDensity0();
	const int numParticles = (int)model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real invH = static_cast<Real>(1.0) / h;
	
	Real density_error = 0.0;

	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Compute pressure accelerations using the current pressure values.
		// (see Algorithm 3, line 7 in [BK17])
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static) 
		for (int i = 0; i < numParticles; i++)
		{
			computePressureAccel(fluidModelIndex, i, density0, m_simulationData.getPressureRho2Data());
		}

		//////////////////////////////////////////////////////////////////////////
		// Update pressure values
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for reduction(+:density_error) schedule(static) 
		for (int j = 0; j < numParticles; j++)
		{
			Real& p_rho2_i = m_simulationData.getPressureRho2(fluidModelIndex, j);

			if (model->getParticleState(j) != ParticleState::Active) {
				p_rho2_i = 0.0;
				continue;
			}
				
			Real aij_pj = compute_aij_pj(fluidModelIndex, j);
			aij_pj *= h * h;

			//////////////////////////////////////////////////////////////////////////
			// Compute source term: s_i = 1 - rho_adv
			// Note: that due to our multiphase handling, the multiplier rho0
			// is missing here
			//////////////////////////////////////////////////////////////////////////
			const Real& densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, j);

			const Real s_i = static_cast<Real>(1.0) - densityAdv;


			//////////////////////////////////////////////////////////////////////////
			// Update the value p/rho^2 (in [BK17] this is kappa/rho):
			// 
			// alpha_i = -1 / (a_ii * rho_i^2)
			// p_rho2_i = (p_i / rho_i^2)
			// 
			// Therefore, the following lines compute the Jacobi iteration:
			// p_i := p_i + 1/a_ii (source_term_i - a_ij * p_j)
			//////////////////////////////////////////////////////////////////////////
			const Real residuum = min(s_i - aij_pj, static_cast<Real>(0.0));     // r = b - A*p
			//p_rho2_i -= residuum * m_simulationData.getFactor(fluidModelIndex, i);

			p_rho2_i = max(p_rho2_i - 0.5 * (s_i - aij_pj) * m_simulationData.getFactor(fluidModelIndex, j), 0.0);

			//////////////////////////////////////////////////////////////////////////
			// Compute the sum of the density errors
			//////////////////////////////////////////////////////////////////////////
			density_error -= density0 * residuum;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Compute the average density error
	//////////////////////////////////////////////////////////////////////////
	avg_density_err = density_error / numParticles;
}

void TimeStepDFSPHbubbleOp::divergenceSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const Real density0 = model->getDensity0();
	const int numParticles = (int)model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const Real invH = static_cast<Real>(1.0) / h;
	
	Real density_error = 0.0;
	
	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Compute pressure accelerations using the current pressure values.
 		// (see Algorithm 2, line 7 in [BK17])
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			if(model->getParticleState(i) == ParticleState::Disabled){
				continue;
			}
			computePressureAccel(fluidModelIndex, i, density0, m_simulationData.getPressureRho2VData());
		}

		//////////////////////////////////////////////////////////////////////////
		// Update pressure 
		//////////////////////////////////////////////////////////////////////////
		#pragma omp for reduction(+:density_error) schedule(static) 
		for (int i = 0; i < numParticles; i++)
		{
			if(model->getParticleState(i) == ParticleState::Disabled){
				continue;
			}

			Real aij_pj = compute_aij_pj(fluidModelIndex, i);
			aij_pj *= h;

			//////////////////////////////////////////////////////////////////////////
			// Compute source term: s_i = -d rho / dt
			//////////////////////////////////////////////////////////////////////////
			const Real& densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
			const Real s_i = -densityAdv;

			//////////////////////////////////////////////////////////////////////////
			// Update the value p/rho^2:
			// 
			// alpha_i = -1 / (a_ii * rho_i^2)
			// pv_rho2_i = (pv_i / rho_i^2)
			// 
			// Therefore, the following line computes the Jacobi iteration:
			// pv_i := pv_i + 1/a_ii (source_term_i - a_ij * pv_j)
			//////////////////////////////////////////////////////////////////////////
			Real& pv_rho2_i = m_simulationData.getPressureRho2_V(fluidModelIndex, i);
			Real residuum = min(s_i - aij_pj, static_cast<Real>(0.0));     // r = b - A*p

			unsigned int numNeighbors = 0;
			// only consider neighbors of same phase here
			numNeighbors += sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i);

			// in case of particle deficiency do not perform a divergence solve
			if (!sim->is2DSimulation())
			{
				if (numNeighbors < 20)
					residuum = 0.0;
			}
			else
			{
				if (numNeighbors < 7)
					residuum = 0.0;
			}
			//pv_rho2_i -= residuum * m_simulationData.getFactor(fluidModelIndex, i);
			pv_rho2_i = max(pv_rho2_i - 0.5*(s_i - aij_pj) * m_simulationData.getFactor(fluidModelIndex, i), 0.0);


			//////////////////////////////////////////////////////////////////////////
			// Compute the sum of the divergence errors
			//////////////////////////////////////////////////////////////////////////
			density_error -= density0 * residuum;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Compute the average divergence error
	//////////////////////////////////////////////////////////////////////////
	avg_density_err = density_error / numParticles;
}

void TimeStepDFSPHbubbleOp::reset()
{
	TimeStep::reset();
	m_nextEmitTime = 0.0;
	m_simulationData.reset();
	m_counter = 0;
	m_iterations = 0;
	m_iterationsV = 0;
	m_numberEmittedTrappedAirParticles = 0;
}

void TimeStepDFSPHbubbleOp::performNeighborhoodSearch()
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

void TimeStepDFSPHbubbleOp::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepDFSPHbubbleOp::resize()
{
	m_simulationData.init();
}

//#ifdef USE_AVX
//
//void TimeStepDFSPH::computeDFSPHFactor(const unsigned int fluidModelIndex)
//{
//	//////////////////////////////////////////////////////////////////////////
//	// Init parameters
//	//////////////////////////////////////////////////////////////////////////
//
//	Simulation* sim = Simulation::getCurrent();
//	const unsigned int nFluids = sim->numberOfFluidModels();
//	FluidModel* model = sim->getFluidModel(fluidModelIndex);
//	const int numParticles = (int)model->numActiveParticles();
//	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
//
//	#pragma omp parallel default(shared)
//	{
//		//////////////////////////////////////////////////////////////////////////
//		// Compute pressure stiffness denominator
//		//////////////////////////////////////////////////////////////////////////
//
//		#pragma omp for schedule(static)
//		for (int i = 0; i < numParticles; i++)
//		{
//			//////////////////////////////////////////////////////////////////////////
//			// Compute gradient dp_i/dx_j * (1/kappa)  and dp_j/dx_j * (1/kappa)
//			// (see Equation (8) and the previous one [BK17])
//			// Note: That in all quantities rho0 is missing due to our
//			// implementation of multiphase simulations.
//			//////////////////////////////////////////////////////////////////////////
//			const Vector3r& xi = model->getPosition(i);
//
//			Real sum_grad_p_k;
//			Vector3r grad_p_i;
//			Vector3f8 xi_avx(xi);
//			Scalarf8 sum_grad_p_k_avx(0.0f);
//			Vector3f8 grad_p_i_avx;
//			grad_p_i_avx.setZero();
//
//			//////////////////////////////////////////////////////////////////////////
//			// Fluid
//			//////////////////////////////////////////////////////////////////////////
//			forall_fluid_neighbors_avx_nox(
//				compute_xj(fm_neighbor, pid);
//				compute_Vj(fm_neighbor);
//				compute_Vj_gradW();
//				sum_grad_p_k_avx += V_gradW.squaredNorm();
//				grad_p_i_avx += V_gradW;
//			);
//
//			//////////////////////////////////////////////////////////////////////////
//			// Boundary
//			//////////////////////////////////////////////////////////////////////////
//			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
//			{
//				forall_boundary_neighbors_avx(
//					const Scalarf8 V_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVolume(0), count);
//					const Vector3f8 grad_p_j = CubicKernel_AVX::gradW(xj_avx - xi_avx) * V_avx;
//					grad_p_i_avx -= grad_p_j;
//				);
//			}
//
//			sum_grad_p_k = sum_grad_p_k_avx.reduce();
//			grad_p_i[0] = grad_p_i_avx.x().reduce();
//			grad_p_i[1] = grad_p_i_avx.y().reduce();
//			grad_p_i[2] = grad_p_i_avx.z().reduce();
//
//			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
//			{
//				forall_density_maps(
//					grad_p_i -= gradRho;
//				);
//			}
//			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
//			{
//				forall_volume_maps(
//					const Vector3r grad_p_j = -Vj * sim->gradW(xi - xj);
//					grad_p_i -= grad_p_j;
//				);
//			}
//
//			sum_grad_p_k += grad_p_i.squaredNorm();
//
//			//////////////////////////////////////////////////////////////////////////
//			// Compute factor alpha_i / rho_i (see Equation (11) in [BK17])
//			//////////////////////////////////////////////////////////////////////////
//			Real& factor = m_simulationData.getFactor(fluidModelIndex, i);
//			if (sum_grad_p_k > m_eps)
//				factor = static_cast<Real>(1.0) / (sum_grad_p_k);
//			else
//				factor = 0.0;
//		}
//	}
//}
//
///** Compute rho_adv,i^(0) (see equation in Section 3.3 in [BK17])
//  * using the velocities after the non-pressure forces were applied.
//**/
//void TimeStepDFSPH::computeDensityAdv(const unsigned int fluidModelIndex, const unsigned int i, const Real h, const Real density0)
//{
//	Simulation *sim = Simulation::getCurrent();
//	FluidModel *model = sim->getFluidModel(fluidModelIndex);
//	const Real &density = model->getDensity(i);
//	Real &densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
//	const Vector3r &xi = model->getPosition(i);
//	const Vector3r &vi = model->getVelocity(i);
//	Real delta = 0.0;
//	const unsigned int nFluids = sim->numberOfFluidModels();
//	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
//
//	Scalarf8 delta_avx(0.0f);
//	const Vector3f8 xi_avx(xi);
//	Vector3f8 vi_avx(vi);
//
//	//////////////////////////////////////////////////////////////////////////
//	// Fluid
//	//////////////////////////////////////////////////////////////////////////
//	forall_fluid_neighbors_avx_nox(
//		compute_xj(fm_neighbor, pid);
//		compute_Vj(fm_neighbor);
//		compute_Vj_gradW();
//		const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &fm_neighbor->getVelocity(0), count);
//		delta_avx += (vi_avx - vj_avx).dot(V_gradW);
//	);
//
//	//////////////////////////////////////////////////////////////////////////
//	// Boundary
//	//////////////////////////////////////////////////////////////////////////
//	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
//	{
//		forall_boundary_neighbors_avx(
//			const Scalarf8 Vj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVolume(0), count);
//			const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVelocity(0), count);
//			delta_avx += Vj_avx * (vi_avx - vj_avx).dot(CubicKernel_AVX::gradW(xi_avx - xj_avx));
//		);
//	}
//
//	delta = delta_avx.reduce();
//
//	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
//	{
//		forall_density_maps(
//			Vector3r vj;
//			bm_neighbor->getPointVelocity(xi, vj);
//			delta -= (vi - vj).dot(gradRho);
//		);
//	}
//	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
//	{
//		forall_volume_maps(
//			Vector3r vj;
//			bm_neighbor->getPointVelocity(xj, vj);
//			delta += Vj * (vi - vj).dot(sim->gradW(xi - xj));
//		);
//	}
//
//	densityAdv = density / density0 + h*delta;
//}
//
///** Compute rho_adv,i^(0) (see equation (9) in Section 3.2 [BK17])
//  * using the velocities after the non-pressure forces were applied.
//  */
//void TimeStepDFSPH::computeDensityChange(const unsigned int fluidModelIndex, const unsigned int i, const Real h)
//{
//	Simulation *sim = Simulation::getCurrent();
//	FluidModel *model = sim->getFluidModel(fluidModelIndex);
//	const Vector3r &xi = model->getPosition(i);
//	const Vector3r &vi = model->getVelocity(i);
//	unsigned int numNeighbors = 0;
//	const unsigned int nFluids = sim->numberOfFluidModels();
//	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
//
//	Scalarf8 densityAdv_avx(0.0f);
//	const Vector3f8 xi_avx(xi);
//	Vector3f8 vi_avx(vi);
//
//	//////////////////////////////////////////////////////////////////////////
//	// Fluid
//	//////////////////////////////////////////////////////////////////////////
//	forall_fluid_neighbors_avx_nox(
//		compute_xj(fm_neighbor, pid);
//		compute_Vj(fm_neighbor);
//		compute_Vj_gradW();
//		const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &fm_neighbor->getVelocity(0), count);
//		densityAdv_avx += (vi_avx - vj_avx).dot(V_gradW);
//	);
//
//	//////////////////////////////////////////////////////////////////////////
//	// Boundary
//	//////////////////////////////////////////////////////////////////////////
//	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
//	{
//		forall_boundary_neighbors_avx(
//			const Scalarf8 Vj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVolume(0), count);
//			const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVelocity(0), count);
//			densityAdv_avx += Vj_avx * (vi_avx - vj_avx).dot(CubicKernel_AVX::gradW(xi_avx - xj_avx));
//		);
//	}
//
//	Real &densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
//	densityAdv = densityAdv_avx.reduce();
//
//	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
//	{
//		forall_density_maps(
//			Vector3r vj;
//			bm_neighbor->getPointVelocity(xi, vj);
//			densityAdv -= (vi - vj).dot(gradRho);
//		);
//	}
//	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
//	{
//		forall_volume_maps(
//			Vector3r vj;
//			bm_neighbor->getPointVelocity(xj, vj);
//			densityAdv += Vj * (vi - vj).dot(sim->gradW(xi - xj));
//		);
//	}
//}
//
///** Compute pressure accelerations using the current pressure values of the particles
// */
//void TimeStepDFSPH::computePressureAccel(const unsigned int fluidModelIndex, const unsigned int i, const Real density0, std::vector<std::vector<Real>>& pressure_rho2, const bool applyBoundaryForces)
//{
//	Simulation* sim = Simulation::getCurrent();
//	FluidModel* model = sim->getFluidModel(fluidModelIndex);
//	const unsigned int nFluids = sim->numberOfFluidModels();
//	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
//
//	Vector3r& ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
//
//	if (model->getParticleState(i) != ParticleState::Active)
//		return;
//
//	// p_rho2_i = (p_i / rho_i^2)
//	const Real p_rho2_i = pressure_rho2[fluidModelIndex][i];
//	const Vector3r &xi = model->getPosition(i);
//
//	Scalarf8 p_rho2_i_avx(p_rho2_i);
//	const Vector3f8 xi_avx(xi);
//	Vector3f8 delta_ai_avx;
//	delta_ai_avx.setZero();
//
//	//////////////////////////////////////////////////////////////////////////
//	// Fluid
//	//////////////////////////////////////////////////////////////////////////
//	forall_fluid_neighbors_avx_nox(
//		compute_xj(fm_neighbor, pid);
//		compute_Vj(fm_neighbor);
//		compute_Vj_gradW();
//		const Scalarf8 densityFrac_avx(fm_neighbor->getDensity0() / density0);
//
//		// p_rho2_j = (p_j / rho_j^2)
//		const Scalarf8 p_rho2_j_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &pressure_rho2[pid][0], count);
//		const Scalarf8 pSum = p_rho2_i_avx + densityFrac_avx * p_rho2_j_avx;
//		delta_ai_avx -= V_gradW * pSum;
//	)
//
//	//////////////////////////////////////////////////////////////////////////
//	// Boundary
//	//////////////////////////////////////////////////////////////////////////
//	if (fabs(p_rho2_i) > m_eps)
//	{
//		if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
//		{
//			const Scalarf8 mi_avx(model->getMass(i));
//			forall_boundary_neighbors_avx(
//				const Scalarf8 Vj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVolume(0), count);
//
//				// Directly update velocities instead of storing pressure accelerations
//				const Vector3f8 a = -CubicKernel_AVX::gradW(xi_avx - xj_avx) * (Vj_avx * p_rho2_i_avx);
//				delta_ai_avx += a;
//
//				if (applyBoundaryForces)
//					bm_neighbor->addForce(xj_avx, -a * mi_avx, count);
//			);
//		}
//	}
//
//	ai[0] = delta_ai_avx.x().reduce();
//	ai[1] = delta_ai_avx.y().reduce();
//	ai[2] = delta_ai_avx.z().reduce();
//
//	if (fabs(p_rho2_i) > m_eps)
//	{
//		if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
//		{
//			forall_density_maps(
//				const Vector3r a = (Real) 1.0 * p_rho2_i * gradRho;
//				ai += a;
//
//				if (applyBoundaryForces)
//					bm_neighbor->addForce(xj, -model->getMass(i) * a);
//			);
//		}
//		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
//		{
//			forall_volume_maps(
//				const Vector3r grad_p_j = -Vj * sim->gradW(xi - xj);
//				const Vector3r a = (Real) 1.0 * p_rho2_i * grad_p_j;
//				ai += a;
//
//				if (applyBoundaryForces)
//					bm_neighbor->addForce(xj, -model->getMass(i) * a);
//			);
//		}
//	}
//}
//
//
//Real TimeStepDFSPH::compute_aij_pj(const unsigned int fluidModelIndex, const unsigned int i)
//{
//	Simulation* sim = Simulation::getCurrent();
//	FluidModel* model = sim->getFluidModel(fluidModelIndex);
//	const unsigned int nFluids = sim->numberOfFluidModels();
//	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
//
//	//////////////////////////////////////////////////////////////////////////
//	// Compute A*p which is the change of the density when applying the
//	// pressure forces.
//	// \sum_j a_ij * p_j = h^2 \sum_j V_j (a_i - a_j) * gradW_ij
//	// This is the RHS of Equation (12) in [BK17]
//	//////////////////////////////////////////////////////////////////////////
//	const Vector3r& xi = model->getPosition(i);
//	const Vector3r& ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
//	const Vector3f8 xi_avx(xi);
//	const Vector3f8 ai_avx(ai);
//	Scalarf8 aij_pj_avx;
//	aij_pj_avx.setZero();
//
//	//////////////////////////////////////////////////////////////////////////
//	// Fluid
//	//////////////////////////////////////////////////////////////////////////
//	forall_fluid_neighbors_avx_nox(
//		compute_xj(fm_neighbor, pid);
//		compute_Vj(fm_neighbor);
//		compute_Vj_gradW();
//
//		const Vector3f8 aj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &m_simulationData.getPressureAccel(pid, 0), count);
//		aij_pj_avx += (ai_avx - aj_avx).dot(V_gradW);
//	);
//
//	//////////////////////////////////////////////////////////////////////////
//	// Boundary
//	//////////////////////////////////////////////////////////////////////////
//	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
//	{
//		forall_boundary_neighbors_avx(
//			const Scalarf8 Vj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVolume(0), count);
//			aij_pj_avx += Vj_avx * ai_avx.dot(CubicKernel_AVX::gradW(xi_avx - xj_avx));
//		);
//	}
//
//	Real aij_pj = aij_pj_avx.reduce();
//
//	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
//	{
//		forall_density_maps(
//			aij_pj -= ai.dot(gradRho);
//		);
//	}
//	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
//	{
//		forall_volume_maps(
//			aij_pj += Vj * ai.dot(sim->gradW(xi - xj));
//		);
//	}
//	return aij_pj;
//}




void TimeStepDFSPHbubbleOp::computeDFSPHFactor(const unsigned int fluidModelIndex)
{
	//////////////////////////////////////////////////////////////////////////
	// Init parameters
	//////////////////////////////////////////////////////////////////////////

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	const int numParticles = (int) model->numActiveParticles();

	#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Compute pressure stiffness denominator
		//////////////////////////////////////////////////////////////////////////

		#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			if (model->getParticleState(i) == ParticleState::Disabled){
				continue;
			}

			//////////////////////////////////////////////////////////////////////////
			// Compute gradient dp_i/dx_j * (1/kappa)  and dp_j/dx_j * (1/kappa)
			// (see Equation (8) and the previous one [BK17])
			// Note: That in all quantities rho0 is missing due to our
			// implementation of multiphase simulations.
			//////////////////////////////////////////////////////////////////////////
			const Vector3r &xi = model->getPosition(i);
			Real sum_grad_p_k = 0.0;
			Vector3r grad_p_i;
			grad_p_i.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Vector3r grad_p_j = -model->getVolume(neighborIndex) * sim->gradW(xi - xj);
				sum_grad_p_k += grad_p_j.squaredNorm();
				grad_p_i -= grad_p_j;
			);
			
			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					const Vector3r grad_p_j = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
					grad_p_i -= grad_p_j;
				);
			}
			
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					grad_p_i -= gradRho;
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					const Vector3r grad_p_j = -Vj * sim->gradW(xi - xj);
					grad_p_i -= grad_p_j;
				);
			}		

			sum_grad_p_k += grad_p_i.squaredNorm();

			//////////////////////////////////////////////////////////////////////////
			// Compute factor as: factor_i = -1 / (a_ii * rho_i^2)
			// where a_ii is the diagonal entry of the linear system 
			// for the pressure A * p = source term
			//////////////////////////////////////////////////////////////////////////
			Real &factor = m_simulationData.getFactor(fluidModelIndex, i);
			if (sum_grad_p_k > m_eps)
				factor = static_cast<Real>(1.0) / (sum_grad_p_k);
			else
				factor = 0.0;
		}
	}
}

/** Compute rho_adv,i^(0) (see equation in Section 3.3 in [BK17])
  * using the velocities after the non-pressure forces were applied.
**/
void TimeStepDFSPHbubbleOp::computeDensityAdv(const unsigned int fluidModelIndex, const unsigned int i, const Real h, const Real density0)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const Real& density = model->getDensity(i);
	Real& densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
	const Vector3r& xi = model->getPosition(i);
	const Vector3r& vi = model->getVelocity(i);
	Real delta = 0.0;
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors_in_same_phase(
		const Vector3r & vj = model->getVelocity(neighborIndex);
		delta += (vi - vj).dot(sim->gradW(xi - xj));
	//delta += fm_neighbor->getVolume(neighborIndex) * (vi - vj).dot(sim->gradW(xi - xj));
	);
	// assumes that all fluid particles have the same volume
	delta *= model->getVolume(i);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors(
			const Vector3r & vj = bm_neighbor->getVelocity(neighborIndex);
		delta += bm_neighbor->getVolume(neighborIndex) * (vi - vj).dot(sim->gradW(xi - xj));
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	{
		forall_density_maps(
			Vector3r vj;
		bm_neighbor->getPointVelocity(xi, vj);
		delta -= (vi - vj).dot(gradRho);
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	{
		forall_volume_maps(
			Vector3r vj;
		bm_neighbor->getPointVelocity(xj, vj);
		delta += Vj * (vi - vj).dot(sim->gradW(xi - xj));
		);
	}

	densityAdv = density / density0 + h*delta;
	}

/** Compute rho_adv,i^(0) (see equation (9) in Section 3.2 [BK17])
  * using the velocities after the non-pressure forces were applied.
  */
void TimeStepDFSPHbubbleOp::computeDensityChange(const unsigned int fluidModelIndex, const unsigned int i, const Real h)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = sim->getFluidModel(fluidModelIndex);
	Real &densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
	const Vector3r &xi = model->getPosition(i);
	const Vector3r& vi = model->getVelocity(i);
	densityAdv = 0.0;
	unsigned int numNeighbors = 0;
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors_in_same_phase(
		const Vector3r & vj = model->getVelocity(neighborIndex);
		densityAdv += (vi - vj).dot(sim->gradW(xi - xj));
	);
	// assumes that all fluid particles have the same volume
	densityAdv *= model->getVolume(i);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors(
			const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
			densityAdv += bm_neighbor->getVolume(neighborIndex) * (vi - vj).dot(sim->gradW(xi - xj));
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	{
		forall_density_maps(
			Vector3r vj;
			bm_neighbor->getPointVelocity(xi, vj);
			densityAdv -= (vi - vj).dot(gradRho);
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	{
		forall_volume_maps(
			Vector3r vj;
			bm_neighbor->getPointVelocity(xj, vj);
			densityAdv += Vj * (vi - vj).dot(sim->gradW(xi - xj));
		);
	}
}

/** Compute pressure accelerations using the current pressure values of the particles
 */
void TimeStepDFSPHbubbleOp::computePressureAccel(const unsigned int fluidModelIndex, const unsigned int i, const Real density0, std::vector<std::vector<Real>>& pressure_rho2, const bool applyBoundaryForces)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	Vector3r& ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
	ai.setZero();

	if (model->getParticleState(i) != ParticleState::Active)
		return;

	// p_rho2_i = (p_i / rho_i^2)
	const Real p_rho2_i = pressure_rho2[fluidModelIndex][i];
	const Vector3r &xi = model->getPosition(i);

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors_in_same_phase(
		// p_rho2_j = (p_j / rho_j^2)
		const Real p_rho2_j = pressure_rho2[fluidModelIndex][neighborIndex];
		const Real pSum = p_rho2_i + p_rho2_j;
		if (fabs(pSum) > m_eps)
		{
			const Vector3r grad_p_j = -model->getVolume(neighborIndex) * sim->gradW(xi - xj);
			ai += pSum * grad_p_j;		
		}
	)

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (fabs(p_rho2_i) > m_eps)
	{
		if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
		{
			forall_boundary_neighbors(
				const Vector3r grad_p_j = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);

				const Vector3r a = (Real) 1.0 * p_rho2_i * grad_p_j;		
				ai += a;
				if (applyBoundaryForces)
					bm_neighbor->addForce(xj, -model->getMass(i) * a);
			);
		}
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		{
			forall_density_maps(
				const Vector3r a = (Real) 1.0 * p_rho2_i * gradRho;			
				ai += a;
				if (applyBoundaryForces)
					bm_neighbor->addForce(xj, -model->getMass(i) * a);
			);
		}
		else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		{
			forall_volume_maps(
				const Vector3r grad_p_j = -Vj * sim->gradW(xi - xj);
				const Vector3r a = (Real) 1.0 * p_rho2_i * grad_p_j;		
				ai += a;

				if (applyBoundaryForces)
					bm_neighbor->addForce(xj, -model->getMass(i) * a);  
			);
		}
	}
}


Real TimeStepDFSPHbubbleOp::compute_aij_pj(const unsigned int fluidModelIndex, const unsigned int i)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	//////////////////////////////////////////////////////////////////////////
	// Compute A*p which is the change of the density when applying the 
	// pressure forces. 
	// \sum_j a_ij * p_j = h^2 \sum_j V_j (a_i - a_j) * gradW_ij
	// This is the RHS of Equation (12) in [BK17]
	//////////////////////////////////////////////////////////////////////////
	const Vector3r& xi = model->getPosition(i);
	const Vector3r& ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
	Real aij_pj = 0.0;

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors_in_same_phase(
		const Vector3r & aj = m_simulationData.getPressureAccel(fluidModelIndex, neighborIndex);
		//aij_pj += fm_neighbor->getVolume(neighborIndex) * (ai - aj).dot(sim->gradW(xi - xj));
		aij_pj += (ai - aj).dot(sim->gradW(xi - xj));
	);
	// assumes that all fluid particles have the same volume
	aij_pj *= model->getVolume(i);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors(
			aij_pj += bm_neighbor->getVolume(neighborIndex) * ai.dot(sim->gradW(xi - xj));
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	{
		forall_density_maps(
			aij_pj -= ai.dot(gradRho);
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	{
		forall_volume_maps(
			aij_pj += Vj * ai.dot(sim->gradW(xi - xj));
		);
	}
	return aij_pj;
}

// Two-way coupling
void TimeStepDFSPHbubbleOp::computeDragForce(const unsigned int fluidModelIndex, const Real h){
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
		Real dragConstant = m_dragConstantAir;
		if (model->getId() != "Air"){
			dragConstant = m_dragConstantLiq;
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

void TimeStepDFSPHbubbleOp::computeBouyancyForce(const unsigned int fluidModelIndex, const Real h){
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

void TimeStepDFSPHbubbleOp::computeCohesionForce(const unsigned int fluidModelIndex, const Real h) {
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

void TimeStepDFSPHbubbleOp::computeViscosityForce(const unsigned int fluidModelIndex, const unsigned int i, const Real h){
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

void TimeStepDFSPHbubbleOp::computeSurfaceTensionForce(const unsigned int fluidModelIndex, const Real h){
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	unsigned int nFluids = sim->numberOfFluidModels();
	unsigned int numParticles = model->numActiveParticles();

	if(model->getId() == "Air") {
		LOG_ERR << "Surface tension force not implemented for air particles.";
		return;
	}

#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
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
}

// #endif

// emits one particle with the given position and velocity direction
void TimeStepDFSPHbubbleOp::emitAirParticleFromVelocityField(unsigned int &numEmittedParticles, const Vector3r vel, const Vector3r pos)//(std::vector <unsigned int> &reusedParticles, unsigned int &indexReuse, unsigned int &numEmittedParticles, const Vector3r& vel, const Vector3r& pos)
{
	// Explanation: TODO
	TimeManager *tm = TimeManager::getCurrent();
	const Real t = tm->getTime();
	const Real timeStepSize = tm->getTimeStepSize();
	Simulation *sim = Simulation::getCurrent();
	const Real radius = sim->getParticleRadius();
	const Real diam = static_cast<Real>(2.0)*radius;

	FluidModel* liquidModel = sim->getFluidModel(0)->getId() == "Liquid" ? sim->getFluidModel(0) : sim->getFluidModel(1);
	FluidModel* airModel = sim->getFluidModel(1)->getId() == "Air" ? sim->getFluidModel(1) : sim->getFluidModel(0);


	if (((airModel->numActiveParticles() + numEmittedParticles) < airModel->numParticles())) // || (reusedParticles.size() > 0))
	{
		const unsigned int indexNotReuse = airModel->numActiveParticles() + numEmittedParticles;

		airModel->getPosition(indexNotReuse) = pos;
		airModel->getVelocity(indexNotReuse) = vel;
		airModel->setParticleState(indexNotReuse, ParticleState::Active);
		airModel->setObjectId(indexNotReuse, 0); //?

		// rng lifetime
		static std::default_random_engine rng;
		std::uniform_real_distribution<Real> lifetime(timeStepSize + m_eps, 1.0);

		Real lifetime_i = lifetime(rng);
		m_simulationData.setLifetime(indexNotReuse, lifetime_i);

		numEmittedParticles++;
		m_numberEmittedTrappedAirParticles++;
	}
}
