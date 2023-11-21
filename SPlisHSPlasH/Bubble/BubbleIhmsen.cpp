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
		// TODO
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


