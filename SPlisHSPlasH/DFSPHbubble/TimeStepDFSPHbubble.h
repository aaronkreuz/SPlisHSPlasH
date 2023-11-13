#ifndef __TimeStepDFSPHbubble_h__
#define __TimeStepDFSPHbubble_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataDFSPHbubble.h"
#include "SPlisHSPlasH/SPHKernels.h"

#define USE_WARMSTART
#define USE_WARMSTART_V

namespace SPH
{
	class SimulationDataDFSPHjs;

	/** \brief This class implements the Divergence-free Smoothed Particle Hydrodynamics approach introduced
	* by Bender and Koschier [BK15,BK17,KBST19].
	*
	* References:
	* - [BK15] Jan Bender and Dan Koschier. Divergence-free smoothed particle hydrodynamics. In ACM SIGGRAPH / Eurographics Symposium on Computer Animation, SCA '15, 147-155. New York, NY, USA, 2015. ACM. URL: http://doi.acm.org/10.1145/2786784.2786796
	* - [BK17] Jan Bender and Dan Koschier. Divergence-free SPH for incompressible and viscous fluids. IEEE Transactions on Visualization and Computer Graphics, 23(3):1193-1206, 2017. URL: http://dx.doi.org/10.1109/TVCG.2016.2578335
	* - [KBST19] Dan Koschier, Jan Bender, Barbara Solenthaler, and Matthias Teschner. Smoothed particle hydrodynamics for physically-based simulation of fluids and solids. In Eurographics 2019 - Tutorials. Eurographics Association, 2019. URL: https://interactivecomputergraphics.github.io/SPH-Tutorial
	*/
	class TimeStepDFSPHbubble : public TimeStep
	{
	protected:
		SimulationDataDFSPHbubble m_simulationData;
		unsigned int m_counter;
		const Real m_eps = static_cast<Real>(1.0e-5);
		bool m_enableDivergenceSolver;
		unsigned int m_iterationsV;
		Real m_maxErrorV;
		unsigned int m_maxIterationsV;

		// BUBBLE related constants
		const Real m_dragConstantAir = static_cast<Real>(8.0);
		const Real m_dragConstatnLiq = static_cast<Real>(3.0);
		const Real m_speedSound = static_cast<Real>(343.0);
		const Real m_kmax = static_cast<Real>(6.0);
		const Real m_minBouyancy = static_cast<Real>(14.0);

		//////////////////////////////////////////////////////////////////////////
		// Forces
		//////////////////////////////////////////////////////////////////////////
		// Liquid forces
		void computeViscosityForce(const unsigned int fluidModelIndex, const unsigned int index, const Real h);

		// Air force
		void computeBouyancyForce(const unsigned int fluidModelIndex, const Real h);
		// -

		void computeDragForce(const unsigned int fluidModelIndex, const Real h);



		void computeDFSPHFactor(const unsigned int fluidModelIndex);
		void pressureSolve();
		void computeKappa(const unsigned int fluidModelIndex, const unsigned int index, const Real h);
		void pressureSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void divergenceSolve();
		void divergenceSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err);
		void computeDensityAdv(const unsigned int fluidModelIndex, const unsigned int index, const Real h, const Real density0);
		void computeDensityChange(const unsigned int fluidModelIndex, const unsigned int index, const Real h);

		void computePressureAccel(const unsigned int fluidModelIndex, const unsigned int i, const Real density0, const std::vector<std::vector<Real>>& pressureList, const bool applyBoundaryForces = false);
		void compute_aii(const unsigned int fluidModelIndex, const unsigned int i, const Real h);
		Real compute_aij_pj(const unsigned int fluidModelIndex, const unsigned int i);
		void computeConstantDensitySourceTerm(const unsigned int fluidModelIndex, const unsigned int i, const Real h);
		void computeDivergenceSourceTerm(const unsigned int fluidModelIndex, const unsigned int i, const Real h);

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

		TimeStepDFSPHbubble();
		virtual ~TimeStepDFSPHbubble(void);

		/** perform a simulation step */
		virtual void step();
		virtual void reset();

		virtual void resize();
	};
}

#endif
