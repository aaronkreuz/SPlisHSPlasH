#ifndef __SimulationDataDFSPHbubbleOp_h__
#define __SimulationDataDFSPHbubbleOp_h__

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
	class SimulationDataDFSPHbubbleOp
	{
		public:
			SimulationDataDFSPHbubbleOp();
			virtual ~SimulationDataDFSPHbubbleOp();

		protected:	
			/** \brief factor \f$\alpha_i\f$ */
			std::vector<std::vector<Real>> m_factor;
			/** \brief advected density */
			std::vector<std::vector<Real>> m_density_adv;

			/** \brief stores \f$\frac{p}{\rho^2}\f$ value of the constant density solver */
			std::vector<std::vector<Real>> m_pressure_rho2;
			/** \brief stores \f$\frac{p}{\rho^2}\f$ value of the divergence solver */
			std::vector<std::vector<Real>> m_pressure_rho2_V;
			std::vector<std::vector<Vector3r>> m_pressureAccel;

			// solely for air particles
			std::vector<unsigned int> m_airParticlesForReuse;
			std::vector<unsigned int> m_onSurfaceAir;
			std::vector<Real> m_lifetimeAir;

			// trappedAir
			std::vector<Vector3r> m_posAirParticles;
			std::vector<Vector3r> m_velAirParticles;


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

			std::vector<std::vector<Real>>& getPressureRho2Data() { return m_pressure_rho2; }
			std::vector<std::vector<Real>>& getPressureRho2VData() { return m_pressure_rho2_V; }
			std::vector<unsigned int>& getParticlesForReuse() {return m_airParticlesForReuse; }

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

			FORCE_INLINE const Real getPressureRho2(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_pressure_rho2[fluidIndex][i];
			}

			FORCE_INLINE Real& getPressureRho2(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_pressure_rho2[fluidIndex][i];
			}

			FORCE_INLINE void setPressureRho2(const unsigned int fluidIndex, const unsigned int i, const Real p)
			{
				m_pressure_rho2[fluidIndex][i] = p;
			}

			FORCE_INLINE const Real getPressureRho2_V(const unsigned int fluidIndex, const unsigned int i) const
			{
				return m_pressure_rho2_V[fluidIndex][i];
			}

			FORCE_INLINE Real& getPressureRho2_V(const unsigned int fluidIndex, const unsigned int i)
			{
				return m_pressure_rho2_V[fluidIndex][i];
			}

			FORCE_INLINE void setPressureRho2_V(const unsigned int fluidIndex, const unsigned int i, const Real p)
			{
				m_pressure_rho2_V[fluidIndex][i] = p;
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

			FORCE_INLINE unsigned int getOnSurface(const unsigned int i) const
			{
				return m_onSurfaceAir[i];
			}

			FORCE_INLINE unsigned int& getOnSurface(const unsigned int i)
			{
				return m_onSurfaceAir[i];
			}

			FORCE_INLINE void setOnSurface(const unsigned int i, const unsigned int val)
			{
				m_onSurfaceAir[i] = val;
			}

			FORCE_INLINE Real getLifetime(const unsigned int i) const
			{
				return m_lifetimeAir[i];
			}

			FORCE_INLINE Real& getLifetime(const unsigned int i)
			{
				return m_lifetimeAir[i];
			}

			FORCE_INLINE void setLifetime(const unsigned int i, const Real val)
			{
				m_lifetimeAir[i] = val;
			}

			FORCE_INLINE std::vector<Vector3r>& getPosAirParticles()
			{
                return m_posAirParticles;
            }

			FORCE_INLINE std::vector<Vector3r> getPosAirParticles() const
			{
                return m_posAirParticles;
            }

			FORCE_INLINE std::vector<Vector3r>& getVelAirParticles()
			{
                return m_velAirParticles;
            }

			FORCE_INLINE std::vector<Vector3r> getVelAirParticles() const
			{
                return m_velAirParticles;
            }
	};
}

#endif