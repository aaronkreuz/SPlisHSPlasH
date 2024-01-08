#include "SPlisHSPlasH/Bubble/BubbleBase.h"

using namespace SPH;
using namespace GenParam;

int BubbleBase::COHESION_COEFFICIENT = -1;
int BubbleBase::DRAG_CONSTANT_AIR = -1;
int BubbleBase::DRAG_CONSTANT_LIQ = -1;
int BubbleBase::BUOYANCY_MAX = -1;
int BubbleBase::BUOYANCY_MIN = -1;

BubbleBase::BubbleBase(FluidModel *model):NonPressureForceBase(model)
{
	m_cohesionCoefficient = static_cast<Real>(12.0); // originally 12.0
	m_dragConstantAir = static_cast<Real>(8.0); // originally 8.0
	m_dragConstantLiq = static_cast<Real>(3.0); // originally 3.0
	m_kmaxBuoyancy = static_cast<Real>(6.0); // originally 6.0
	m_minBuoyancy = static_cast<Real>(1.4); // originally 14.0
}

BubbleBase::~BubbleBase(void)
{
}

void BubbleBase::initParameters(){
	NonPressureForceBase::initParameters();

	if(m_model->getId() == "Air"){
		COHESION_COEFFICIENT = createNumericParameter("cohesionCoefficient", "Cohesion coefficient", &m_cohesionCoefficient);
		setGroup(COHESION_COEFFICIENT, "Fluid Model|Bubble Forces");
		setDescription(COHESION_COEFFICIENT, "Cohesion coefficient");

		DRAG_CONSTANT_AIR = createNumericParameter("dragConstantAir", "Drag constant air", &m_dragConstantAir);
		setGroup(DRAG_CONSTANT_AIR, "Fluid Model|Bubble Forces");
		setDescription(DRAG_CONSTANT_AIR, "Drag constant air");

		BUOYANCY_MAX = createNumericParameter("buoyancyMax", "Buoyancy max", &m_kmaxBuoyancy);
		setGroup(BUOYANCY_MAX, "Fluid Model|Bubble Forces");
		setDescription(BUOYANCY_MAX, "Buoyancy max");

		BUOYANCY_MIN = createNumericParameter("buoyancyMin", "Buoyancy min", &m_minBuoyancy);
		setGroup(BUOYANCY_MIN, "Fluid Model|Bubble Forces");
		setDescription(BUOYANCY_MIN, "Buoyancy min");
	}

	if(m_model->getId() == "Liquid"){
		DRAG_CONSTANT_LIQ = createNumericParameter("dragConstantLiq", "Drag constant liquid", &m_dragConstantLiq);
		setGroup(DRAG_CONSTANT_LIQ, "Fluid Model|Bubble Forces");
		setDescription(DRAG_CONSTANT_LIQ, "Drag constant liquid");
	}
}