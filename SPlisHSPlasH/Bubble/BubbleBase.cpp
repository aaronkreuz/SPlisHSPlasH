#include "BubbleBase.h"

using namespace SPH;
using namespace GenParam;

int BubbleBase::COHESION_COEFFICIENT = -1;

BubbleBase::BubbleBase(FluidModel *model) :
	NonPressureForceBase(model)
{
	m_cohesionCoefficient = static_cast<Real>(12.0);
}

BubbleBase::~BubbleBase(void)
{
}

void BubbleBase::initParameters(void){
	NonPressureForceBase::initParameters();

	COHESION_COEFFICIENT = createNumericParameter("cohesion", "Cohesion coefficient", &m_cohesionCoefficient);
	setGroup(COHESION_COEFFICIENT, "Fluid Model|Cohesion force");
	setDescription(COHESION_COEFFICIENT, "Coefficient for the cohesion force computation");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(COHESION_COEFFICIENT));
	rparam->setMinValue(0.0);
}