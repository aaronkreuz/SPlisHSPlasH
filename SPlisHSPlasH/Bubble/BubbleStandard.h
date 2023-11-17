#ifndef __BubbleStandard_h__
#define __BubbleStandard_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/Bubble/BubbleBase.h"

namespace SPH
{
	/** \brief This class implements the bubble cohesion method introduced
	* by Ihmsen et al. \cite Ihmsen:2014.
	*/
	class BubbleStandard : public BubbleBase
	{
	protected:
		void computeCohesionStandard();
		void computeBouyancyStandard();
		void computeDragOnLiquidStandard();
		void computeDragOnAirStandard();

	public:
		BubbleStandard(FluidModel *model);
		virtual ~BubbleStandard(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new BubbleStandard(model); }

		virtual void step();
		virtual void reset();

		//virtual void performNeighborhoodSearchSort();

	};
}

#endif
