#ifndef __BubbleBase_h__
#define __BubbleBase_h__

#include "SPlisHSPlasH/Common.h"
#include "../FluidModel.h"
#include "../NonPressureForceBase.h"

namespace SPH
{
	/** \brief Base class for all bubble methods.
	*/
	class BubbleBase : public NonPressureForceBase
	{
	protected:
		Real m_cohesionCoefficient;
		Real m_dragConstantAir;
		Real m_dragConstantLiq;
		const Real m_eps = static_cast<Real>(1.0e-5);
		const Real m_speedSound = static_cast<Real>(343.0);
		Real m_kmaxBuoyancy;
		Real m_minBuoyancy;

		virtual void initParameters();

	public:
		static int COHESION_COEFFICIENT;
		static int DRAG_CONSTANT_AIR;
		static int DRAG_CONSTANT_LIQ;
		static int BUOYANCY_MAX;
		static int BUOYANCY_MIN;


		BubbleBase(FluidModel* model);
		virtual ~BubbleBase(void);
	};
}

#endif