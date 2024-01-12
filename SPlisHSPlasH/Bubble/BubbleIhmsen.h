#ifndef __BubbleIhmsen_h__
#define __BubbleIhmsen_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/Bubble/BubbleBase.h"

namespace SPH
{
	/** \brief This class implements provides non-pressure force functions for the implementation of the bubble-framework
	* by Ihmsen et al. \cite Ihmsen:2014.
	*
	* Currently it is important that there are two fluids, one for the liquid and one for the air. The fluid IDs are currently hardcoded.
	* "Air" and "Liquid".
	*/
	class BubbleIhmsen : public BubbleBase
	{
	protected:
		std::vector<Vector3r> m_normals;
		std::vector<Matrix3r> m_L_air; // kernel gradient correction matrix for air particles
		// std::vector<unsigned int> m_onSurface; // Moved to simulationData
		// std::vector<Real> m_lifetime; // Moved to simulationData

		Real m_onSurfaceThresholdDensity = 0.5; // TODO: changeable

		bool m_useGradientCorrection = false;
		int m_cohesionForce = 0;
		int m_buoyancyForce = 0;
		int m_dragForceLiq = 0;
		int m_dragForceAir = 0;

		void computeForces(FluidModel* model);

		void computeCohesionIhmsen(FluidModel* model);
		void computeCohesionIhmsenKernel(FluidModel* model);
		void computeCohesionAkinci2013(FluidModel* model);
		void computeGradientCorrection(void);
		void computeNormals(void);

		void computeBuoyancyIhmsen(FluidModel* model);

		void computeDragIhmsen(FluidModel* model);

		virtual void initParameters();

	public:
		BubbleIhmsen(FluidModel *model);
		virtual ~BubbleIhmsen(void);

		static NonPressureForceBase* creator(FluidModel* model) { return new BubbleIhmsen(model); }

		virtual void step();
		virtual void reset();

		static int COHESION_FORCE;
		static int ENUM_COHESION_FORCE_NONE;
		static int ENUM_COHESION_FORCE_IHMSEN;
		static int ENUM_COHESION_FORCE_SURFACE_TENSION;
		static int ENUM_COHESION_FORCE_IHMSEN_KERNEL;
		static int ENUM_COHESION_FORCE_AKINCI2013;
		static int ENUM_USE_SURFACE_TENSION;

		static int BUOYANCY_FORCE;
		static int ENUM_BUOYANCY_FORCE_NONE;
		static int ENUM_BUOYANCY_FORCE_IHMSEN;

		static int DRAG_FORCE_LIQ;
		static int ENUM_DRAG_FORCE_LIQ_NONE;
		static int ENUM_DRAG_FORCE_LIQ_IHMSEN;

		static int DRAG_FORCE_AIR;
		static int ENUM_DRAG_FORCE_AIR_NONE;
		static int ENUM_DRAG_FORCE_AIR_IHMSEN;

		static int USE_GRADIENT_CORRECTION;

		FORCE_INLINE Vector3r &getNormal(const unsigned int i)
		{
			return m_normals[i];
		}
		
		FORCE_INLINE const Vector3r &getNormal(const unsigned int i) const
		{
			return m_normals[i];
		}
		
		FORCE_INLINE void setNormal(const unsigned int i, const Vector3r &val)
		{
			m_normals[i] = val;
		}

		virtual void performNeighborhoodSearchSort();

	};
}

#endif
