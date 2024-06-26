#ifndef __SimulationDataPCISPHtest_h__
#define __SimulationDataPCISPHtest_h__

#include "SPlisHSPlasH/Common.h"
#include <vector>
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH
{

	class SimulationDataPCISPHtest
	{
	public:
		SimulationDataPCISPHtest();
		virtual ~SimulationDataPCISPHtest();

	protected:
		// INFO: mass and density are stored in the fluid model

		// INFO: Accessing all arrays via m_array[fluidModelIndex][particleIndex]

		std::vector<std::vector<Vector3r>> m_predX;
		std::vector<std::vector<Vector3r>> m_predV;

		// density advection (predicted density)
		std::vector<std::vector<Real>> m_densityAdv; 

		std::vector<std::vector<Real>> m_pressure;
		std::vector<std::vector<Vector3r>> m_pressureAccel;

		// scaling factor for pressure
		std::vector<Real> m_pcisph_factor;

	public:
		/** Initialize the arrays containing the particle data.
		*/
		virtual void init();

		/** Release the arrays containing the particle data.
		*/
		virtual void cleanup();

		/** Reset the particle data.
		*/
		virtual void reset();

		/** Important: First call m_model->performNeighborhoodSearchSort()
		 * to call the z_sort of the neighborhood search.
		 */
		void performNeighborhoodSearchSort();

		void emittedParticles(FluidModel* model, const unsigned int startIndex);

		Real getPcisphFactor(const unsigned int fluidIndex) {
            return m_pcisph_factor[fluidIndex];
        }

		FORCE_INLINE Vector3r& getPredictedPosition(const unsigned int fluidIndex, const unsigned int i)
		{
            return m_predX[fluidIndex][i];
        }

		FORCE_INLINE const Vector3r& getPredictedPosition(const unsigned int fluidIndex, const unsigned int i) const
		{
            return m_predX[fluidIndex][i];
        }

		FORCE_INLINE void setPredictedPosition(const unsigned int fluidIndex, const unsigned int i, const Vector3r& pos)
		{
            m_predX[fluidIndex][i] = pos;
        }

		FORCE_INLINE Vector3r& getPredictedVelocity(const unsigned int fluidIndex, const unsigned int i) {
            return m_predV[fluidIndex][i];
        }

		FORCE_INLINE const Vector3r& getPredictedVelocity(const unsigned int fluidIndex, const unsigned int i) const
		{
            return m_predV[fluidIndex][i];
        }

		FORCE_INLINE void setPredictedVelocity(const unsigned int fluidIndex, const unsigned int i, const Vector3r& vel)
		{
            m_predV[fluidIndex][i] = vel;
        }

		FORCE_INLINE Real& getDensityAdv(const unsigned int fluidIndex, const unsigned int i) {
            return m_densityAdv[fluidIndex][i];
        }

		FORCE_INLINE const Real getDensityAdv(const unsigned int fluidIndex, const unsigned int i) const
		{
            return m_densityAdv[fluidIndex][i];
        }

		FORCE_INLINE void setDensityAdv(const unsigned int fluidIndex, const unsigned int i, const Real d)
		{
            m_densityAdv[fluidIndex][i] = d;
        }

		FORCE_INLINE const Real getPressure(const unsigned int fluidIndex, const unsigned int i) const
		{
			return m_pressure[fluidIndex][i];
		}

		FORCE_INLINE Real& getPressure(const unsigned int fluidIndex, const unsigned int i)
		{
			return m_pressure[fluidIndex][i];
		}

		FORCE_INLINE void setPressure(const unsigned int fluidIndex, const unsigned int i, const Real p)
		{
			m_pressure[fluidIndex][i] = p;
		}

		FORCE_INLINE Vector3r& getPressureAccel(const unsigned int fluidIndex, const unsigned int i)
		{
			return m_pressureAccel[fluidIndex][i];
		}

		FORCE_INLINE const Vector3r& getPressureAccel(const unsigned int fluidIndex, const unsigned int i) const
		{
			return m_pressureAccel[fluidIndex][i];
		}

		FORCE_INLINE void setPressureAccel(const unsigned int fluidIndex, const unsigned int i, const Vector3r& val)
		{
			m_pressureAccel[fluidIndex][i] = val;
		}

	};
}

#endif