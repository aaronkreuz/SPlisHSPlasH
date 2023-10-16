#ifndef __TimeStepPCISPHtest_h__
#define __TimeStepPCISPHtest_h__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/TimeStep.h"
#include "SimulationDataPCISPHtest.h"
#include "SPlisHSPlasH/SPHKernels.h"

namespace SPH
{
	class SimulationDataPCISPHtest;

	/** \brief This class implements the Weakly Compressible SPH for Free Surface Flows approach introduced
	* by Becker and Teschner [BT07].
	*
	* References:
	* - [BT07] Markus Becker and Matthias Teschner. Weakly compressible SPH for free surface flows. In ACM SIGGRAPH/Eurographics Symposium on Computer Animation, SCA '07, 209-217. Aire-la-Ville, Switzerland, Switzerland, 2007. Eurographics Association. URL: http://dl.acm.org/citation.cfm?id=1272690.1272719
	*/
	class TimeStepPCISPHtest : public TimeStep
	{
	protected:
		SimulationDataPCISPHtest m_simulationData;
		unsigned int m_counter;


		/** Perform the neighborhood search for all fluid particles.
		*/
		void performNeighborhoodSearch();

		virtual void emittedParticles(FluidModel* model, const unsigned int startIndex);
		virtual void initParameters();

	public:
		TimeStepPCISPHtest();
		virtual ~TimeStepPCISPHtest(void);

		virtual void step();
		virtual void reset();
		virtual void resize();
	};
}

#endif
