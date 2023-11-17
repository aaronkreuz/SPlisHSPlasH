#include "BubbleStandard.h"
#include <iostream>
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/Utilities/MathFunctions.h"

using namespace SPH;

BubbleStandard::BubbleStandard(FluidModel *model) :
	BubbleBase(model)
{
    m_cohesionCoefficient = static_cast<Real>(12.0);
}

BubbleStandard::~BubbleStandard(void)
{
}

void BubbleStandard::step()
{
    FluidModel *model = m_model;

	std::cout << "Perform step" << std::endl;


	//////////////////////////////////////////////////////////////////////////
	// TODO
	//////////////////////////////////////////////////////////////////////////
}

void BubbleStandard::reset()
{
}

void computeCohesionStandard(){

}

void computeBouyancyStandard(){

}

void computeDragOnLiquidStandard(){

}

void computeDragOnAirStandard(){

}


