#include "BubbleIhmsen.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"

using namespace SPH;
using namespace GenParam;

int BubbleIhmsen::COHESION_FORCE = -1;
int BubbleIhmsen::ENUM_COHESION_FORCE_NONE = -1;
int BubbleIhmsen::ENUM_COHESION_FORCE_IHMSEN = -1;
int BubbleIhmsen::ENUM_COHESION_FORCE_SURFACE_TENSION = -1;
int BubbleIhmsen::ENUM_COHESION_FORCE_IHMSEN_KERNEL = -1;

int BubbleIhmsen::BUOYANCY_FORCE = -1;
int BubbleIhmsen::ENUM_BUOYANCY_FORCE_NONE = -1;
int BubbleIhmsen::ENUM_BUOYANCY_FORCE_IHMSEN = -1;

int BubbleIhmsen::DRAG_FORCE_LIQ = -1;
int BubbleIhmsen::ENUM_DRAG_FORCE_LIQ_NONE = -1;
int BubbleIhmsen::ENUM_DRAG_FORCE_LIQ_IHMSEN = -1;

int BubbleIhmsen::DRAG_FORCE_AIR = -1;
int BubbleIhmsen::ENUM_DRAG_FORCE_AIR_NONE = -1;
int BubbleIhmsen::ENUM_DRAG_FORCE_AIR_IHMSEN = -1;


BubbleIhmsen::BubbleIhmsen(FluidModel *model) :
	BubbleBase(model)
{
	m_normals.resize(model->numParticles(), Vector3r::Zero());
	model->addField({ "normal", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_normals[i][0]; }, false });

	if(model->getId() == "Air"){
		m_cohesionForce = 1;
		m_buoyancyForce = 1;
		m_dragForceLiq = 0;
		m_dragForceAir = 1;
	}
	else if(model->getId() == "Liquid"){
		m_cohesionForce = 0;
		m_buoyancyForce = 0;
		m_dragForceLiq = 1;
		m_dragForceAir = 0;
	}

}

BubbleIhmsen::~BubbleIhmsen(void)
{
	m_model->removeFieldByName("normal");
	m_normals.clear();
}

void BubbleIhmsen::initParameters()
{
	BubbleBase::initParameters();

	if(m_model->getId() == "Air"){
		COHESION_FORCE = createEnumParameter("cohesionForceType", "Cohesion force", &m_cohesionForce);
		setGroup(COHESION_FORCE, "Fluid Model|Bubble Forces");
		setDescription(COHESION_FORCE, "Method for the cohesion force computation.");
		EnumParameter* enumParam = static_cast<EnumParameter*>(getParameter(COHESION_FORCE));
		enumParam->addEnumValue("None", ENUM_COHESION_FORCE_NONE);
		enumParam->addEnumValue("Ihmsen", ENUM_COHESION_FORCE_IHMSEN);
		enumParam->addEnumValue("Surface Tension", ENUM_COHESION_FORCE_SURFACE_TENSION);
		enumParam->addEnumValue("Ihmsen Kernel", ENUM_COHESION_FORCE_IHMSEN_KERNEL);

		BUOYANCY_FORCE = createEnumParameter("buoyancyForceType", "Buoyancy force", &m_buoyancyForce);
		setGroup(BUOYANCY_FORCE, "Fluid Model|Bubble Forces");
		setDescription(BUOYANCY_FORCE, "buoyancy force computation.");
		enumParam = static_cast<EnumParameter*>(getParameter(BUOYANCY_FORCE));
		enumParam->addEnumValue("None", ENUM_BUOYANCY_FORCE_NONE);
		enumParam->addEnumValue("Ihmsen", ENUM_BUOYANCY_FORCE_IHMSEN);

		DRAG_FORCE_AIR = createEnumParameter("dragForceTypeAir", "Drag force air", &m_dragForceAir);
		setGroup(DRAG_FORCE_AIR, "Fluid Model|Bubble Forces");
		setDescription(DRAG_FORCE_AIR, "Method for the drag force computation on air particles.");
		enumParam = static_cast<EnumParameter*>(getParameter(DRAG_FORCE_AIR));
		enumParam->addEnumValue("None", ENUM_DRAG_FORCE_AIR_NONE);
		enumParam->addEnumValue("Ihmsen", ENUM_DRAG_FORCE_AIR_IHMSEN);
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

	//////////////////////////////////////////////////////////////////////////
	// Buoyancy
	//////////////////////////////////////////////////////////////////////////
	if(m_buoyancyForce == ENUM_BUOYANCY_FORCE_IHMSEN)
	{
		computeBouyancyIhmsen(model);
	}

	//////////////////////////////////////////////////////////////////////////
	// Drag
	//////////////////////////////////////////////////////////////////////////
	if(m_dragForceAir == ENUM_DRAG_FORCE_AIR_IHMSEN)
	{
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
			const Real densj = model->getDensity(neighborIndex);
			const Vector3r xij = xi - xj;
			acc_cohesion += densj * xij * sim->W(xij);
			);

		acceleration -= m_cohesionCoefficient * acc_cohesion;
	 }
}

void BubbleIhmsen::computeCohesionAkinci2013(FluidModel* model){
	// CohesionKernel::W();
	Simulation* sim = Simulation::getCurrent();
	unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const Real supportRadius = sim->getSupportRadius();
	const unsigned int numParticles = model->numActiveParticles();
	const Real density0 = m_model->getDensity0();

	computeNormals();

	for(int i = 0; i < numParticles; i++){
		Vector3r& accel = model->getAcceleration(i);
		Vector3r a_ij = Vector3r::Zero();
		const Vector3r& xi = model->getPosition(i);
		const Vector3r& ni = getNormal(i);
		const Real& rhoi = m_model->getDensity(i);

		forall_fluid_neighbors_in_same_phase(
			Vector3r xij = (xi-xj);
			const Real mj = model->getMass(neighborIndex);
			const Real& rhoj = model->getDensity(neighborIndex);
			const Real K_ij = static_cast<Real>(2.0)*density0 / (rhoi + rhoj);
			a_ij.setZero();

			a_ij -= m_cohesionCoefficient * mj * (xij / xij.norm()) * CohesionKernel::W(xi-xj);

			const Vector3r& nj = getNormal(neighborIndex);
			a_ij -= m_cohesionCoefficient * (ni-nj);

			accel += K_ij * a_ij;
		);
	}

	// Should I consider the boundary here?

}

// Standard method used in BUBBLE-paper
void BubbleIhmsen::computeBouyancyIhmsen(FluidModel* model){
	 Simulation* sim = Simulation::getCurrent();
	 const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	 unsigned int numParticles = model->numActiveParticles();
	 const Vector3r grav(sim->getVecValue<Real>(Simulation::GRAVITATION));

	 for(int i = 0; i < numParticles; i++){
		Vector3r& acceleration = model->getAcceleration(i);

		Vector3r acc_bouyancy = Vector3r::Zero();

		// look for number of air particle neighbors
		int numNeighbors = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i);

		// equation (10)
		acc_bouyancy = m_minBuoyancy * (m_kmaxBuoyancy - (m_kmaxBuoyancy - 1) * exp(-0.1 * numNeighbors)) * grav;
		acceleration -= acc_bouyancy;
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

			const Real pi_ij = std::max(0.0f, ((vi - vj).dot(xi - xj))/((xi - xj).norm() + m_eps * (h2)));
			// const Real pi_ij = max(0.0f, ((vi - vj).dot(xi - xj))/((xi - xj).squaredNorm() + m_eps * (h2)));

			acc_drag += m_j * ((dragConstant*h*m_speedSound)/(density_j + density_i)) * pi_ij * sim->gradW(xi - xj);
			);

		acceleration += acc_drag;
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
			const Vector3r &xi = m_model->getPosition(i);
			Vector3r &ni = getNormal(i);
			ni.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Real density_j = m_model->getDensity(neighborIndex);
				ni += m_model->getMass(neighborIndex) / density_j * sim->gradW(xi - xj);
				)
				ni = supportRadius*ni;
		}
	 }

}
