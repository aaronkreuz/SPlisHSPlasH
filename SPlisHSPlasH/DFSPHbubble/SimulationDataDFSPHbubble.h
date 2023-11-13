#ifndef __SimulationDataDFSPHbubble_h__
#define __SimulationDataDFSPHbubble_h__

#include "SPlisHSPlasH/Common.h"
#include <vector>
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH 
{	
	/** \brief Simulation data which is required by the method Divergence-free Smoothed Particle Hydrodynamics introduced
	* by Bender and Koschier [BK15,BK17].
	*
	* References:
	* - [BK15] Jan Bender and Dan Koschier. Divergence-free smoothed particle hydrodynamics. In ACM SIGGRAPH / Eurographics Symposium on Computer Animation, SCA '15, 147-155. New York, NY, USA, 2015. ACM. URL: http://doi.acm.org/10.1145/2786784.2786796
	* - [BK17] Jan Bender and Dan Koschier. Divergence-free SPH for incompressible and viscous fluids. IEEE Transactions on Visualization and Computer Graphics, 23(3):1193-1206, 2017. URL: http://dx.doi.org/10.1109/TVCG.2016.2578335
	*/
	class SimulationDataDFSPHbubble
	{
		public:
			SimulationDataDFSPHbubble();
			virtual ~SimulationDataDFSPHbubble();

		protected:	
			/** \brief factor \f$\alpha_i\f$ */
			std::vector<std::vector<Real>> m_factor;
			/** \brief advected density */
			std::vector<std::vector<Real>> m_density_adv;

			std::vector<std::vector<Real>> m_aii;

			// constant density source term
			std::vector<std::vector<Real>> m_source_term;

			// divergence source term
			std::vector<std::vector<Real>> m_source_term_div;

			std::vector<std::vector<Real>> m_pressure;
			std::vector<std::vector<Real>> m_pressure_V;

			std::vector<std::vector<Vector3r>> m_pressureAccel;

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
			void emittedParticles(FluidModel *model, const unsigned int startIndex);

			std::vector<std::vector<Real>>& getPressureRho2Data() { return m_pressure; }
			std::vector<std::vector<Real>>& getPressureRho2VData() { return m_pressure_V; }

			FORCE_INLINE const Real getSourceTermDiv(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_source_term_div[fluidIndex][i];
			}

			FORCE_INLINE Real& getSourceTermDiv(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_source_term_div[fluidIndex][i];
			}

			FORCE_INLINE void setSourceTermDiv(const unsigned int fluidIndex, const unsigned int i, const Real s)
			{
				m_source_term_div[fluidIndex][i] = s;
			}

			FORCE_INLINE const Real getSourceTerm(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_source_term[fluidIndex][i];
			}

			FORCE_INLINE Real& getSourceTerm(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_source_term[fluidIndex][i];
			}

			FORCE_INLINE void setSourceTerm(const unsigned int fluidIndex, const unsigned int i, const Real s){
				m_source_term[fluidIndex][i] = s;
			}

			FORCE_INLINE const Real getDiagElement(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_aii[fluidIndex][i];
			}

			FORCE_INLINE Real& getDiagElement(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_aii[fluidIndex][i];
			}

			FORCE_INLINE void setDiagElement(const unsigned int fluidIndex, const unsigned int i, const Real aii)
			{
				m_aii[fluidIndex][i] = aii;
			}

			FORCE_INLINE std::vector<Real>& getFactorsModel(const unsigned int fluidIndex) {
				return m_factor[fluidIndex];
			}

			FORCE_INLINE const Real getFactor(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_factor[fluidIndex][i];
			}

			FORCE_INLINE Real& getFactor(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_factor[fluidIndex][i];
			}

			FORCE_INLINE void setFactor(const unsigned int fluidIndex, const unsigned int i, const Real p)
			{
				m_factor[fluidIndex][i] = p;
			}

			FORCE_INLINE const Real getDensityAdv(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_density_adv[fluidIndex][i];
			}

			FORCE_INLINE Real& getDensityAdv(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_density_adv[fluidIndex][i];
			}

			FORCE_INLINE void setDensityAdv(const unsigned int fluidIndex, const unsigned int i, const Real d)
			{
				m_density_adv[fluidIndex][i] = d;
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

			FORCE_INLINE const Real getPressure_V(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_pressure_V[fluidIndex][i];
			}

			FORCE_INLINE Real& getPressure_V(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_pressure_V[fluidIndex][i];
			}

			FORCE_INLINE void setPressure_V(const unsigned int fluidIndex, const unsigned int i, const Real p)
			{
				m_pressure_V[fluidIndex][i] = p;
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