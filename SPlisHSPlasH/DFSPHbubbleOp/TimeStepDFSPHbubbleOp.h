#ifndef __TimeStepDFSPHbubbleOp_h__
#define __TimeStepDFSPHbubbleOp_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataDFSPHbubbleOp.h"
#include "SPlisHSPlasH/SPHKernels.h"

#define USE_WARMSTART
#define USE_WARMSTART_V

namespace SPH
{
	class SimulationDataDFSPH;

	/** \brief This class implements the Divergence-free Smoothed Particle Hydrodynamics approach introduced
	* by Bender and Koschier [BK15,BK17,KBST19].
	*
	* References:
	* - [BK15] Jan Bender and Dan Koschier. Divergence-free smoothed particle hydrodynamics. In ACM SIGGRAPH / Eurographics Symposium on Computer Animation, SCA '15, 147-155. New York, NY, USA, 2015. ACM. URL: http://doi.acm.org/10.1145/2786784.2786796
	* - [BK17] Jan Bender and Dan Koschier. Divergence-free SPH for incompressible and viscous fluids. IEEE Transactions on Visualization and Computer Graphics, 23(3):1193-1206, 2017. URL: http://dx.doi.org/10.1109/TVCG.2016.2578335
	* - [KBST19] Dan Koschier, Jan Bender, Barbara Solenthaler, and Matthias Teschner. Smoothed particle hydrodynamics for physically-based simulation of fluids and solids. In Eurographics 2019 - Tutorials. Eurographics Association, 2019. URL: https://interactivecomputergraphics.github.io/SPH-Tutorial
	*/
	class TimeStepDFSPHbubbleOp : public TimeStep
	{
	protected:
		SimulationDataDFSPHbubbleOp m_simulationData;
		unsigned int m_counter;
		const Real m_eps = static_cast<Real>(1.0e-5);
		bool m_enableDivergenceSolver;
		unsigned int m_iterationsV;
		Real m_maxErrorV;
		unsigned int m_maxIterationsV;
		bool m_enableTrappedAir;
		bool m_enableTrappedAirOptimization; // do not emit air particle if another air particle is too close (threshold)
		bool m_enableIsolationCriterion; // reduce lifetime of particles if they have less than threshold neighbors

		unsigned int m_iterationsLiq;
		unsigned int m_iterationsAir;
		unsigned int m_iterationsVliq;
		unsigned int m_iterationsVair;

		// BUBBLE related constants
		Real m_dragConstantAir = static_cast<Real>(8.0);
		Real m_dragConstantLiq = static_cast<Real>(3.0);
		const Real m_speedSound = static_cast<Real>(343.0);
		Real m_kmax = static_cast<Real>(6.0);
		Real m_minBouyancy = static_cast<Real>(2.4);
		const Real m_cohesionConstant = static_cast<Real>(12.0);
		const Real m_surfaceTensionConstant = static_cast<Real>(1.0);
		Real m_onSurfaceThresholdDensity = static_cast<Real>(0.3);

		// Trapped air
		int m_trappedAirApproach = 0; // 0: Ihmsen et al. 2011, 1: Ihmsen et al. 2012
		Real m_initialEmitTime = static_cast<Real>(0.1); // AK 2024: avoid emitting in the first timestep
		Real m_nextEmitTime = static_cast<Real>(0.1); // AK 2024: avoid emitting in the first timestep
		Real m_emitTimeDistance = static_cast<Real>(0.1); // AK 2024
		unsigned int m_maxAirParticlesPerTimestep = static_cast<unsigned int>(20);

		Real m_vMinTrappedAir = 9.0f; // Ihmsen et al. 2011 and 2012
		Real m_vtTrappedAir = 0.3f; // Ihmsen et al. 2011

		Real m_vDiffThresholdMin = static_cast<Real>(5.0); // Ihmsen et al. 2012
		Real m_vDiffThresholdMax = static_cast<Real>(20.0); // Ihmsen et al. 2012

		Real m_thresholdCavitationDensityRatio = static_cast<Real>(0.5); // AK 2024

		int m_numberEmittedTrappedAirParticles = 0;

		// different maxError in pressure solver for Air phase
		Real m_maxErrorAir; // AK 2024


		//////////////////////////////////////////////////////////////////////////
		// Forces
		//////////////////////////////////////////////////////////////////////////
		// Liquid forces
		void computeViscosityForce(const unsigned int fluidModelIndex, const unsigned int index, const Real h);
		void computeSurfaceTensionForce(const unsigned int fluidModelIndex, const Real h);

		// Air force
		void computeBouyancyForce(const unsigned int fluidModelIndex, const Real h);

		void computeCohesionForce(const unsigned int fluidModelIndex, const Real h);
		// -

		void computeDragForce(const unsigned int fluidModelIndex, const Real h);

		void computeDFSPHFactor(const unsigned int fluidModelIndex);
		void pressureSolve();
		void pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void divergenceSolve();
		void divergenceSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void computeDensityAdv(const unsigned int fluidModelIndex, const unsigned int index, const Real h, const Real density0);
		void computeDensityChange(const unsigned int fluidModelIndex, const unsigned int index, const Real h);

		void computePressureAccel(const unsigned int fluidModelIndex, const unsigned int i, const Real density0, std::vector<std::vector<Real>>& pressure_rho2, const bool applyBoundaryForces = false);
		Real compute_aij_pj(const unsigned int fluidModelIndex, const unsigned int i);

		void computeOnSurfaceAir();
		void trappedAirIhmsen2011(const unsigned int fluidModelIndex, const unsigned int i, unsigned int &numTrappedAirParticles, std::vector<unsigned int>& indicesGen);
		void trappedAirIhmsen2012(unsigned int& emittedParticles, std::vector<unsigned int>& indicesGen);
		void trappedAirCavitation(const unsigned int fluidModelIndex, const unsigned int i, unsigned int& numTrappedAirParticles, std::vector<unsigned int>& indicesGen);

		void emitAirParticleFromVelocityField(unsigned int &numEmittedParticles, Vector3r vel, Vector3r pos);

		/** Perform the neighborhood search for all fluid particles.
		*/
		void performNeighborhoodSearch();
		virtual void emittedParticles(FluidModel *model, const unsigned int startIndex);

		/** Init all generic parameters */
		virtual void initParameters();

	public:
		static int SOLVER_ITERATIONS_V;
		static int MAX_ITERATIONS_V;
		static int MAX_ERROR_V;
		static int USE_DIVERGENCE_SOLVER;
		static int SOLVER_ITERATIONS_LIQ;
		static int SOLVER_ITERATIONS_AIR;
		static int SOLVER_ITERATIONS_V_LIQ;
		static int SOLVER_ITERATIONS_V_AIR;
		static int MAX_ERROR_AIR;
		static int SURFACE_THRESHOLD_DENSITY_RATIO;
		static int ENABLE_ISOLATION_CRITERION;

		// trapped air
		static int TRAPPED_AIR_APPROACH;
		static int ENUM_TRAPPED_AIR_APPROACH_NONE;
		static int ENUM_TRAPPED_AIR_APPROACH_IHMSEN2011;
		static int ENUM_TRAPPED_AIR_APPROACH_IHMSEN2012;
		static int ENUM_TRAPPED_AIR_APPROACH_CAVITATION;

		static int USE_TRAPPED_AIR;
		static int USE_TRAPPED_AIR_OPTIMIZATION;
		static int VMIN_TRAPPED_AIR;
		static int VT_TRAPPED_AIR;
		static int VDIFF_THRESHOLD_MIN; // Ihmsen et al. 2012
		static int VDIFF_THRESHOLD_MAX; // Ihmsen et al. 2012
		static int MAX_AIR_PARTICLES_PER_STEP;
		static int EMIT_TIME_DISTANCE;
		static int DENSITY_RATIO_CAVITATION; // AK 2024
		static int NEXT_EMIT_TIME; // AK 2024
		static int INITIAL_EMIT_TIME; // AK 2024



		TimeStepDFSPHbubbleOp();
		virtual ~TimeStepDFSPHbubbleOp(void);

		FORCE_INLINE unsigned int getOnSurface(const unsigned int i) const {
			return m_simulationData.getOnSurface(i);
		}
		FORCE_INLINE void setOnSurface(const unsigned int i, const unsigned int val) {
			m_simulationData.setOnSurface(i, val);
		}

		/** perform a simulation step */
		virtual void step();
		virtual void reset();

		virtual void resize();
	};
}

#endif
