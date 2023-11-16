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

		// BUBBLE related constants
		Real m_dragConstantAir = static_cast<Real>(8.0);
		Real m_dragConstantLiq = static_cast<Real>(3.0);
		const Real m_speedSound = static_cast<Real>(343.0);
		Real m_kmax = static_cast<Real>(6.0);
		Real m_minBouyancy = static_cast<Real>(2.4);
		const Real m_cohesionConstant = static_cast<Real>(12.0);
		const Real m_surfaceTensionConstant = static_cast<Real>(1.0);

		int m_cohesionForce;

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

		static int DRAG_COEFFICIENT_AIR;
		static int DRAG_COEFFICIENT_LIQ;

		static int COHESION_FORCE_TYPE;
		static int ENUM_COHESION_FORCE_NONE;
		static int ENUM_COHESION_FORCE_IHMSEN;
		static int ENUM_COHESION_FORCE_SURFACE_TENSION;

		TimeStepDFSPHbubbleOp();
		virtual ~TimeStepDFSPHbubbleOp(void);

		/** perform a simulation step */
		virtual void step();
		virtual void reset();

		virtual void resize();
	};
}

#endif