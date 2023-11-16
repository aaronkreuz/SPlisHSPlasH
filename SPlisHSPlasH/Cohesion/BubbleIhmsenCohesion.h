#ifndef __BubbleIhmsenCohesion_h__
#define __BubbleIhmsenCohesion_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	/** \brief This class implements the bubble cohesion method introduced
	* by Ihmsen et al. \cite Ihmsen:2014.
	*/
	class BubbleIhmsenCohesion : public NonPressureForceBase
	{
	protected:
        Real m_cohesionCoefficient;

		virtual void initParameters();

	public:
		BubbleIhmsenCohesion(FluidModel *model);
		virtual ~BubbleIhmsenCohesion(void);

        static NonPressureForceBase* creator(FluidModel* model) { return new BubbleIhmsenCohesion(model); }

		virtual void step();
		virtual void reset();

		virtual void performNeighborhoodSearchSort();

        static int COHESION_COEFFICIENT;

	};
}

#endif
