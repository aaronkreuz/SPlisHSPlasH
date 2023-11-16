#include "BubbleIhmsenCohesion.h"
#include <iostream>
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/Utilities/MathFunctions.h"

using namespace SPH;
using namespace std;

int BubbleIhmsenCohesion::COHESION_COEFFICIENT = -1;

BubbleIhmsenCohesion::BubbleIhmsenCohesion(FluidModel *model) :
	NonPressureForceBase(model)
{
    m_cohesionCoefficient = static_cast<Real>(12.0);
}

BubbleIhmsenCohesion::~BubbleIhmsenCohesion(void)
{
}

void BubbleIhmsenCohesion::step()
{
    FluidModel *model = m_model;
}

void BubbleIhmsenCohesion::reset()
{
}

void BubbleIhmsenCohesion::initParameters()
{
    NonPressureForceBase::initParameters();
	
    COHESION_COEFFICIENT= createNumericParameter("cohesionCoefficient", "Cohesion force coefficient", &m_cohesionCoefficient);
	setGroup(COHESION_COEFFICIENT, "Simulation|BUBBLE");
	setDescription(COHESION_COEFFICIENT, "Coefficient for the cohesion force.");
}


