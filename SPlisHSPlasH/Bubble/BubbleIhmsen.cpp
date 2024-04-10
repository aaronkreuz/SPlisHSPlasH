#include "BubbleIhmsen.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/DFSPHbubbleOp/TimeStepDFSPHbubbleOp.h"
#include "SPlisHSPlasH/Utilities/MathFunctions.h"


using namespace SPH;
using namespace GenParam;

int BubbleIhmsen::COHESION_FORCE = -1;
int BubbleIhmsen::ENUM_COHESION_FORCE_NONE = -1;
int BubbleIhmsen::ENUM_COHESION_FORCE_IHMSEN = -1;
int BubbleIhmsen::ENUM_COHESION_FORCE_SURFACE_TENSION = -1;
int BubbleIhmsen::ENUM_COHESION_FORCE_IHMSEN_KERNEL = -1;
int BubbleIhmsen::ENUM_COHESION_FORCE_AKINCI2013 = -1;
int BubbleIhmsen::ENUM_USE_SURFACE_TENSION = -1;

int BubbleIhmsen::BUOYANCY_FORCE = -1;
int BubbleIhmsen::ENUM_BUOYANCY_FORCE_NONE = -1;
int BubbleIhmsen::ENUM_BUOYANCY_FORCE_IHMSEN = -1;
int BubbleIhmsen::ENUM_BUOYANCY_FORCE_DISPLACEMENT = -1;

int BubbleIhmsen::DRAG_FORCE_LIQ = -1;
int BubbleIhmsen::ENUM_DRAG_FORCE_LIQ_NONE = -1;
int BubbleIhmsen::ENUM_DRAG_FORCE_LIQ_IHMSEN = -1;

int BubbleIhmsen::DRAG_FORCE_AIR = -1;
int BubbleIhmsen::ENUM_DRAG_FORCE_AIR_NONE = -1;
int BubbleIhmsen::ENUM_DRAG_FORCE_AIR_IHMSEN = -1;

int BubbleIhmsen::USE_GRADIENT_CORRECTION = -1;

int BubbleIhmsen::COHESION_COEFFICIENT_NORMAL = -1;


BubbleIhmsen::BubbleIhmsen(FluidModel *model) :
	BubbleBase(model)
{
	m_normals.resize(model->numParticles(), Vector3r::Zero());
	if (model->getId() == "Air") {
		m_L_air.resize(model->numParticles(), Matrix3r::Identity());
    }	
	// model->addField({ "normal", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_normals[i][0]; }, false });

	if(model->getId() == "Air"){
		m_cohesionForce = 5;
		m_buoyancyForce = 1;
		m_dragForceLiq = 0;
		m_dragForceAir = 1;
		m_cohesionCoefficientNormal = 5.0;
	}
	else if(model->getId() == "Liquid"){
		m_cohesionForce = 0;
		m_buoyancyForce = 0;
		m_dragForceLiq = 1;
		m_dragForceAir = 0;
	}

	m_useGradientCorrection = false;
	m_onSurfaceThresholdDensity = 0.5; // TODO: changeable
}

BubbleIhmsen::~BubbleIhmsen(void)
{
	// m_model->removeFieldByName("normal");
	m_normals.clear();
	m_L_air.clear();
}

void BubbleIhmsen::initParameters()
{
	BubbleBase::initParameters();

	if(m_model->getId() == "Air"){
		COHESION_COEFFICIENT_NORMAL = createNumericParameter("cohesionCoefficientNormal", "Cohesion coefficient normal", &m_cohesionCoefficientNormal);
		setGroup(COHESION_COEFFICIENT_NORMAL, "Fluid Model|Bubble Forces");
		setDescription(COHESION_COEFFICIENT_NORMAL, "Cohesion coefficient for the normal part of the cohesion force (only used for Akinci et al. 2013).");

		COHESION_FORCE = createEnumParameter("cohesionForceType", "Cohesion force", &m_cohesionForce);
		setGroup(COHESION_FORCE, "Fluid Model|Bubble Forces");
		setDescription(COHESION_FORCE, "Method for the cohesion force computation.");
		EnumParameter* enumParam = static_cast<EnumParameter*>(getParameter(COHESION_FORCE));
		enumParam->addEnumValue("None", ENUM_COHESION_FORCE_NONE);
		enumParam->addEnumValue("Ihmsen", ENUM_COHESION_FORCE_IHMSEN);
		enumParam->addEnumValue("Surface Tension", ENUM_COHESION_FORCE_SURFACE_TENSION);
		enumParam->addEnumValue("Ihmsen Kernel", ENUM_COHESION_FORCE_IHMSEN_KERNEL);
		enumParam->addEnumValue("Akinci 2013", ENUM_COHESION_FORCE_AKINCI2013);
		enumParam->addEnumValue("Use Surface Tension methods", ENUM_USE_SURFACE_TENSION);

		BUOYANCY_FORCE = createEnumParameter("buoyancyForceType", "Buoyancy force", &m_buoyancyForce);
		setGroup(BUOYANCY_FORCE, "Fluid Model|Bubble Forces");
		setDescription(BUOYANCY_FORCE, "buoyancy force computation.");
		enumParam = static_cast<EnumParameter*>(getParameter(BUOYANCY_FORCE));
		enumParam->addEnumValue("None", ENUM_BUOYANCY_FORCE_NONE);
		enumParam->addEnumValue("Ihmsen", ENUM_BUOYANCY_FORCE_IHMSEN);
		enumParam->addEnumValue("Displacement", ENUM_BUOYANCY_FORCE_DISPLACEMENT);

		DRAG_FORCE_AIR = createEnumParameter("dragForceTypeAir", "Drag force air", &m_dragForceAir);
		setGroup(DRAG_FORCE_AIR, "Fluid Model|Bubble Forces");
		setDescription(DRAG_FORCE_AIR, "Method for the drag force computation on air particles.");
		enumParam = static_cast<EnumParameter*>(getParameter(DRAG_FORCE_AIR));
		enumParam->addEnumValue("None", ENUM_DRAG_FORCE_AIR_NONE);
		enumParam->addEnumValue("Ihmsen", ENUM_DRAG_FORCE_AIR_IHMSEN);

		USE_GRADIENT_CORRECTION = createBoolParameter("useGradientCorrection", "Use gradient correction", &m_useGradientCorrection);
		setGroup(USE_GRADIENT_CORRECTION, "Fluid Model|Bubble Forces");
		setDescription(USE_GRADIENT_CORRECTION, "Use gradient correction for the drag force computation on air particles.");
	}

	else if(m_model->getId() == "Liquid"){
		DRAG_FORCE_LIQ = createEnumParameter("dragForceTypeLiq", "Drag force liquid", &m_dragForceLiq);
		setGroup(DRAG_FORCE_LIQ, "Fluid Model|Bubble Forces");
		setDescription(DRAG_FORCE_LIQ, "Method for the drag force computation on liquid particles.");
		EnumParameter* enumParam = static_cast<EnumParameter*>(getParameter(DRAG_FORCE_LIQ));
		enumParam->addEnumValue("None", ENUM_DRAG_FORCE_LIQ_NONE);
		enumParam->addEnumValue("Ihmsen", ENUM_DRAG_FORCE_LIQ_IHMSEN);
	}
}


void BubbleIhmsen::step()
{
    FluidModel *model = m_model;

	computeForces(model);
}

void BubbleIhmsen::reset()
{

}

void BubbleIhmsen::computeForces(FluidModel* model){
	//////////////////////////////////////////////////////////////////////////
	// Cohesion
	//////////////////////////////////////////////////////////////////////////
	if(m_cohesionForce == ENUM_COHESION_FORCE_IHMSEN)
	{
		computeCohesionIhmsen(model);
	}
	else if(m_cohesionForce == ENUM_COHESION_FORCE_SURFACE_TENSION)
	{
		// TODO
	}
	else if(m_cohesionForce == ENUM_COHESION_FORCE_IHMSEN_KERNEL)
	{
		computeCohesionIhmsenKernel(model);
	}
	else if (m_cohesionForce == ENUM_COHESION_FORCE_AKINCI2013) 
	{
		computeCohesionAkinci2013(model);
	}
	else if(m_cohesionForce == ENUM_USE_SURFACE_TENSION)
	{
		// just links to the registered surface tension method as cohesion method
		model->computeSurfaceTension();
	}

	//////////////////////////////////////////////////////////////////////////
	// Buoyancy
	//////////////////////////////////////////////////////////////////////////
	if(m_buoyancyForce == ENUM_BUOYANCY_FORCE_IHMSEN)
	{
		computeBuoyancyIhmsen(model);
	}
	else if (m_buoyancyForce == ENUM_BUOYANCY_FORCE_DISPLACEMENT)
	{
		computeBuoyancyDisplacement(model);
	}

	//////////////////////////////////////////////////////////////////////////
	// Drag
	//////////////////////////////////////////////////////////////////////////
	if(m_dragForceAir == ENUM_DRAG_FORCE_AIR_IHMSEN)
	{
		if(m_useGradientCorrection){
            computeGradientCorrection();
        }

		computeDragIhmsen(model);
	}

	if(m_dragForceLiq == ENUM_DRAG_FORCE_LIQ_IHMSEN)
	{
		computeDragIhmsen(model);
	}
}

// Standard method used in BUBBLE-paper
void BubbleIhmsen::computeCohesionIhmsen(FluidModel* model){
	 Simulation* sim = Simulation::getCurrent();
	 const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	 unsigned int numParticles = model->numActiveParticles();

	 for (int i = 0; i < numParticles; i++){

		Vector3r& acceleration = model->getAcceleration(i);
		const Vector3r& xi = model->getPosition(i);
		Vector3r acc_cohesion = Vector3r::Zero();

		forall_fluid_neighbors_in_same_phase(
			if(model->getParticleState(neighborIndex) == ParticleState::Disabled){
				continue;
			}

			const Real densj = model->getDensity(neighborIndex);
			acc_cohesion += densj*(xi - xj);
		);

		acceleration -= m_cohesionCoefficient * acc_cohesion;
	 }
}

// Standard method of Ihmsen et al. supplemented with a Kernel scaling
void BubbleIhmsen::computeCohesionIhmsenKernel(FluidModel* model){
	 Simulation* sim = Simulation::getCurrent();
	 const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	 unsigned int numParticles = model->numActiveParticles();
	 for (int i = 0; i < numParticles; i++){
		Vector3r& acceleration = model->getAcceleration(i);
		const Vector3r& xi = model->getPosition(i);
		Vector3r acc_cohesion = Vector3r::Zero();

		forall_fluid_neighbors_in_same_phase(
			if(model->getParticleState(neighborIndex) == ParticleState::Disabled){
				continue;
			}

			const Real densj = model->getDensity(neighborIndex);
			const Vector3r xij = xi - xj;
			acc_cohesion += densj * xij * sim->W(xij);
			);

		acceleration -= m_cohesionCoefficient * acc_cohesion;
	 }
}

// Own Cohesion Kernel variant, constructed according to Liu and Liu 2010
Real BubbleIhmsen::cohesionKernelAdapted(const Vector3r r)
{
	Simulation* sim = Simulation::getCurrent();
	const Real h = sim->getSupportRadius();
	const Real h2 = h*h;
	const Real h4 = h2*h2;

	const Real factor = 32.0 / (M_PI * h4 * h4 * h);

	const Real dist = r.norm();
	const Real dist3 = dist*dist*dist;
	const Real h_r = h - (dist);

	if (2 * dist > h && dist <= h) {
		return factor * (2 * dist3 * (h_r * h_r * h_r));
	}
	else if (dist > 0 && 2 * dist <= h) {
		return factor * ((-1/32)*(0.5*h - dist) * (0.5 * h - dist) * (0.5 * h - dist) + (2 * dist3* (h_r * h_r * h_r)));
	}

	return 0;
}

// See: "Versatile surface tension and adhesion for SPH fluids" Akinci et al. 2013
void BubbleIhmsen::computeCohesionAkinci2013(FluidModel* model){
	Simulation* sim = Simulation::getCurrent();
	unsigned int fluidModelIndex = model->getPointSetIndex();
	const Real supportRadius = sim->getSupportRadius();
	const unsigned int numParticles = model->numActiveParticles();
	const Real density0 = model->getDensity0();

	computeNormals();

	for(int i = 0; i < numParticles; i++){
		Vector3r& accel = model->getAcceleration(i);
		Vector3r a_ij = Vector3r::Zero();
		const Vector3r& xi = model->getPosition(i);
		const Vector3r& ni = getNormal(i);
		const Real& rhoi = m_model->getDensity(i);

		forall_fluid_neighbors_in_same_phase(
			if(model->getParticleState(neighborIndex) == ParticleState::Disabled){
				continue;
			}

			// Cohesion
			Vector3r xij = (xi-xj);
			const Real mj = model->getMass(neighborIndex);
			const Real& rhoj = model->getDensity(neighborIndex);
			const Real K_ij = static_cast<Real>(2.0)*density0 / (rhoi + rhoj);
			a_ij.setZero();

			a_ij -= m_cohesionCoefficient * mj * (xij / xij.norm()) * CohesionKernel::W(xi-xj);

			// Curvature
			const Vector3r& nj = getNormal(neighborIndex);
			a_ij -= m_cohesionCoefficientNormal * (ni-nj);

			accel += K_ij * a_ij;
		);
	}

	// TODO: Should I consider the boundary here?
}

void SPH::BubbleIhmsen::computeGradientCorrection(void)
{
	assert(m_model->getId() == "Air");
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = m_model; // naming convention, air model
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned nParticles = m_model->numActiveParticles();
	const unsigned nFluids = sim->numberOfFluidModels();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int i = 0; i < nParticles; i++) {
			Matrix3r& L = m_L_air[i];
			L.setZero();

			const Vector3r& xi = model->getPosition(i);

			forall_fluid_neighbors(
				const Vector3r xij = xi - xj;
			const Vector3r xji = xj - xi;

			// const Real m_j = model->getMass(neighborIndex);
			// const Real rho_j = model->getDensity(neighborIndex);
			const Real volume_j = fm_neighbor->getVolume(neighborIndex);

			L += volume_j * sim->gradW(xij) * xji.transpose();
			);

			bool invertible = false;
			L.computeInverseWithCheck(m_L_air[i], invertible, static_cast<Real>(1e-9));
			if (!invertible)
			{
				MathFunctions::pseudoInverse(L, m_L_air[i]);
				//m_L[i] = Matrix3r::Identity();
			}
			L = L.inverse(); // Might be inefficient
		}
	}
	
}

// Standard method used in BUBBLE-paper
void BubbleIhmsen::computeBuoyancyIhmsen(FluidModel* model){
	Simulation* sim = Simulation::getCurrent();
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	unsigned int numParticles = model->numActiveParticles();
	const Vector3r grav(sim->getVecValue<Real>(Simulation::GRAVITATION));

	// TODO: Better solution?
	TimeStepDFSPHbubbleOp* timeStep = (TimeStepDFSPHbubbleOp*)sim->getTimeStep();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for(int i = 0; i < numParticles; i++){
			Vector3r& acceleration = model->getAcceleration(i);

			// check if particle in on the surface and if yes: compute different buoyancy
			if(timeStep->getOnSurface(i)){
				// TODO: Decide whether considering garvity here or not
				// NOTE: Ihmsen oiriginal: acceleration -= grav; -> Does not work. Particles will not stick at fluid surface and fly around.
				// acceleration -= 0.5*grav; 
				continue;
			}

			Vector3r acc_bouyancy = Vector3r::Zero();

			// look for number of air particle neighbors
			int numNeighbors = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i);

			// equation (10)
			acc_bouyancy = m_minBuoyancy * (m_kmaxBuoyancy - (m_kmaxBuoyancy - 1) * exp(-0.1 * numNeighbors)) * grav;
			acceleration -= acc_bouyancy;
		}
	}
}

// Experimental 
void SPH::BubbleIhmsen::computeBuoyancyDisplacement(FluidModel* model)
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	unsigned int numParticles = model->numActiveParticles();
	const Vector3r grav(sim->getVecValue<Real>(Simulation::GRAVITATION));

	// TODO: Better solution?
	TimeStepDFSPHbubbleOp* timeStep = (TimeStepDFSPHbubbleOp*)sim->getTimeStep();

	for (int i = 0; i < numParticles; i++) {

		const Real mass_i = model->getMass(i);
		const Real density_i = model->getDensity(i);

		Real volumeEstimate = mass_i / density_i;

		// estimate volume of bubble of particle i
		forall_fluid_neighbors_in_same_phase(
			const Real mass_j = model->getMass(neighborIndex);
			const Real density_j = model->getDensity(neighborIndex);

			volumeEstimate += mass_j / density_j;
		);

		Vector3r& acceleration = model->getAcceleration(i);

		acceleration -= m_kmaxBuoyancy * (1000.0 * volumeEstimate * grav);
	}
}




// Standard method used in BUBBLE-paper
void BubbleIhmsen::computeDragIhmsen(FluidModel* model){
	 Simulation* sim = Simulation::getCurrent();
	 TimeManager* tm = TimeManager::getCurrent();
	 const Real h = tm->getTimeStepSize();
	 const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	 unsigned int nFluids = sim->numberOfFluidModels();
	 unsigned int numParticles = model->numActiveParticles();
	 const Real h2 = h*h;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
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
				if(fm_neighbor->getParticleState(neighborIndex) == ParticleState::Disabled){
					continue;
				}

				const Real m_j = fm_neighbor->getMass(neighborIndex);
				const Real density_j = fm_neighbor->getDensity(neighborIndex);
				const Vector3r& vj = fm_neighbor->getVelocity(neighborIndex);

				// take the speed of sound in the respective fluid into account
				Real speedSoundFluid = m_speedSoundAir;
				if (fm_neighbor->getId() == "Liquid") {
					speedSoundFluid = m_speedSoundWater;
				}

				const Real pi_ij = std::max(0.0f, ((vi - vj).dot(xi - xj))/((xi - xj).norm() + m_eps * (h2)));
				// const Real pi_ij = max(0.0f, ((vi - vj).dot(xi - xj))/((xi - xj).squaredNorm() + m_eps * (h2)));

				Vector3r gradKernel = sim->gradW(xi - xj);
				if (m_useGradientCorrection && model->getId() == "Air") {
					// Kernel gradient correction
                    gradKernel = m_L_air[i] * gradKernel;
                }

				acc_drag += m_j * ((dragConstant*h* speedSoundFluid)/(density_j + density_i)) * pi_ij * gradKernel;
			);

			acceleration += acc_drag;
		}
	}
}

// copied from SurfaceTension_Akinci2013.cpp
void BubbleIhmsen::computeNormals()
{
	 Simulation *sim = Simulation::getCurrent();
	 const Real supportRadius = sim->getSupportRadius();
	 const unsigned int numParticles = m_model->numActiveParticles();
	 const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	 const unsigned int nFluids = sim->numberOfFluidModels();
	 FluidModel *model = m_model;

	 // Compute normals
	#pragma omp parallel default(shared)
	 {
		#pragma omp for schedule(static)
		for (int i = 0; i < (int)numParticles; i++)
		{
			 if(model->getParticleState(i) == ParticleState::Disabled){
				 continue;
			 }

			const Vector3r &xi = m_model->getPosition(i);
			Vector3r &ni = getNormal(i);
			ni.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Real density_j = m_model->getDensity(neighborIndex);
				const Vector3r xij = xi - xj;
				
				ni += (m_model->getMass(neighborIndex) / (density_j)) * sim->gradW(xij);
			);
			
			ni *= supportRadius;
		}
	 }
}

// !! Moved to TimeStepDFSPHbubbleOp !!
// compute state of the air particles: on surface or inside liquid
// void BubbleIhmsen::computeOnSurface(){
// 	Simulation *sim = Simulation::getCurrent();
// 	const unsigned int numParticles = m_model->numActiveParticles();
// 	const unsigned int nFluids = sim->numberOfFluidModels();
// 	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
// 	FluidModel* model = m_model; // air model
// 	const Vector3r grav(sim->getVecValue<Real>(Simulation::GRAVITATION));
// 	TimeManager* tm = TimeManager::getCurrent();
// 	const Real h = tm->getTimeStepSize();
//
// 	for(int i = 0; i < numParticles; i++){
// 		if(model->getParticleState(i) == ParticleState::Disabled){
// 			continue;
// 		}
//
// 		const Real density_i = m_model->getDensity(i);
// 		const Vector3r& xi = m_model->getPosition(i);
//
// 		// This condition seems to be error prone. If an air particle gets "trapped" it might not have any air-neighbors and would be falsly identified as "on the surface"
// 		// if(density_i > m_onSurfaceThresholdDensity){
// 		// 	m_onSurface[i] = 1;
// 		// }
//
// 		// look for number of liquid particle neighbors of the air particle
// 		// int numLiqNeighbors = sim->numberOfNeighbors(m_model->getPointSetIndex(), fluidModelIndex, i);
// 		// std::cout << numLiqNeighbors << std::endl;
//
// 		volatile bool onSurface = true;
// 		// looping over liquid neighbors
// 		forall_fluid_neighbors_in_different_phase(
// 			if(!onSurface){
// 				continue;
// 			}
// 			if((xj-xi).dot(grav) < 0){
// 				onSurface = false;
// 			}
// 		);
//
// 		m_onSurface[i] = onSurface;
//
// 		if(onSurface){
// 			m_lifetime[i] -= h;
//
// 			if(m_lifetime[i] <= 0.0){
// 				// Disable an air particle at the end of its lifetime
// 				model->setParticleState(i, ParticleState::Disabled);
// 				// -> clean-up in TimeStep
// 				m_onSurface[i] = 0;
// 			}
// 		}
//
// 		// all particles of a bubble should disperse at once.
// 		forall_fluid_neighbors_in_same_phase(
// 			m_lifetime[i] = std::min(m_lifetime[i], m_lifetime[neighborIndex]);
// 		);
// 	}
// }

void BubbleIhmsen::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numActiveParticles();
	if (numPart == 0)
		return;

	Simulation* sim = Simulation::getCurrent();
	auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
	d.sort_field(&m_normals[0]); 
}
