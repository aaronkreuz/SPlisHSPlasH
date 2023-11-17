#ifndef __BubbleBase_h__
#define __BubbleBase_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/NonPressureForceBase.h"

namespace SPH
{
	/** \brief Base class for all cohesion methods.
	 *  At the moment solely used for Air-particles.
	*/
	class BubbleBase : public NonPressureForceBase
	{
	protected:
		Real m_cohesionCoefficient;

		virtual void initParameters();

	public:
		static int COHESION_COEFFICIENT;

		BubbleBase(FluidModel* model);
		virtual ~BubbleBase(void);
	};
}

#endif