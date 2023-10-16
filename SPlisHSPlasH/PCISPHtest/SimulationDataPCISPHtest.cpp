#include "SimulationDataPCISPHtest.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "../Simulation.h"

using namespace SPH;

SimulationDataPCISPHtest::SimulationDataPCISPHtest()
{
}

SimulationDataPCISPHtest::~SimulationDataPCISPHtest(void)
{
	cleanup();
}

// INFO: Initialize arrays with correct size
// INFO: Needs to be called at first in TimeStep constructor
void SimulationDataPCISPHtest::init()
{
	// DEBUG
	std::cout << "init simulation data" << std::endl;
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	m_predX.resize(nModels);
	m_predV.resize(nModels);
	m_densityAdv.resize(nModels);

	m_pressure.resize(nModels);
	m_pressureAccel.resize(nModels);

	m_pcisph_factor.resize(nModels);

	// iterate over fluid models and resize arrays according to number of particles in each model
	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel* fm = sim->getFluidModel(i);

		m_predX[i].resize(fm->numParticles(), fm->getPosition0(i)); // TODO: is this correct?
		m_predV[i].resize(fm->numParticles(), Vector3r::Zero()); // TODO: intialize with initial velocity?
		m_densityAdv[i].resize(fm->numParticles(), 0.0f);

		m_pressure[i].resize(fm->numParticles(), 0.0f);
		m_pressureAccel[i].resize(fm->numParticles(), Vector3r::Zero());
	}

	// initialize pcisph constant per particle
	// use prototype particle with max. fluid neighbors (index 0 in this case)
	for (auto i = 0; i < nModels; i++) {
		FluidModel* fm = sim->getFluidModel(i);

		auto rho0 = fm->getDensity0();

		// see Solenthaler and Pajarola [SP09]. Note that squared timestep will be added in TimeStepPCISPHtest::step().
		const Real beta = 2.0f * fm->getVolume(0) * fm->getVolume(0);
		
		// compute sum of gradients of kernel function
		Vector3r sumGradW = Vector3r::Zero();
		
		// compute squared sum of gradients of kernel function
		Real sumGradW2 = 0.0;

		const auto supportRadius = sim->getSupportRadius();
		const auto squaredSupportRadius = supportRadius * supportRadius;
		const auto particleRadius = sim->getParticleRadius();
		const auto diam = static_cast<Real>(2.0) * particleRadius;

		// regular sampling around (0,0,0)
		auto xi = Vector3r(0,0,0);
		// actual computation of gradient sums
		for (auto j = 0; j < fm->numActiveParticles(); j++) {
			
			// if sim is 2D
			if (sim->is2DSimulation()) {
				auto xj = Vector3r(-supportRadius, -supportRadius, 0);

				while (xj[0] <= supportRadius) {
					while (xj[1] <= supportRadius) {
						const auto squaredDist = (xj[0] - xi[0]) * (xj[0] - xi[0]) + (xj[1] - xi[1]) * (xj[1] - xi[1]);

						// check if xj is in the support of xi
						if (squaredDist <= squaredSupportRadius) {

							auto gradW = sim->gradW(xi - xj);
							sumGradW += gradW;
							sumGradW2 += gradW.squaredNorm();
						}
						xj[1] += diam;
					}
					xj[0] += diam;
					xj[1] = -supportRadius;
				}
			}

			// 3D
			else {
				auto xj = Vector3r(-supportRadius, -supportRadius, -supportRadius);
				while(xj[0] <= supportRadius) {
					while(xj[1] <= supportRadius) {
						while(xj[2] <= supportRadius) {
                            const auto squaredDist = (xj[0] - xi[0]) * (xj[0] - xi[0]) + (xj[1] - xi[1]) * (xj[1] - xi[1]) + (xj[2] - xi[2]) * (xj[2] - xi[2]);

                            // check if xj is in the support of xi
							if (squaredDist <= squaredSupportRadius) {

                                auto gradW = sim->gradW(xi - xj);
                                sumGradW += gradW;
                                sumGradW2 += gradW.squaredNorm();
                            }
                            xj[2] += diam;
                        }
                        xj[1] += diam;
                        xj[2] = -supportRadius;
                    }
                    xj[0] += diam;
                    xj[1] = -supportRadius;
                }
			}

		}

		// finally compute pcisph factor
		m_pcisph_factor[i] = static_cast<Real>(1.0) / (beta * (sumGradW.squaredNorm() + sumGradW2));

	}
}

void SimulationDataPCISPHtest::cleanup()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		m_predX[i].clear();
		m_predV[i].clear();
		m_densityAdv[i].clear();

		m_pressure[i].clear();
		m_pressureAccel[i].clear();
	}
	m_predX.clear();
	m_predV.clear();
	m_densityAdv.clear();

	m_pressure.clear();
	m_pressureAccel.clear();
}

// INFO: Reset arrays to initial values
void SimulationDataPCISPHtest::reset()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel* fm = sim->getFluidModel(i);
		for (unsigned int j = 0; j < fm->numActiveParticles(); j++)
		{
			m_predX[i][j] = fm->getPosition0(j); // TODO: is this correct?
			m_predV[i][j].setZero(); // TODO: intialize with initial velocity?
			m_densityAdv[i][j] = 0.0f; // TODO: float?
			m_pressure[i][j] = 0.0f; // TODO: float?
			m_pressureAccel[i][j].setZero();
		}
	}
}

// TODO: What is this for? What is exactly happening here?
void SimulationDataPCISPHtest::performNeighborhoodSearchSort()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel* fm = sim->getFluidModel(i);
		const unsigned int numPart = fm->numActiveParticles();
		if (numPart != 0)
		{
			CompactNSearch::PointSet const& d = sim->getNeighborhoodSearch()->point_set(fm->getPointSetIndex());

			d.sort_field(&m_predX[i][0]);
			d.sort_field(&m_predV[i][0]);
			d.sort_field(&m_densityAdv[i][0]);
			d.sort_field(&m_pressure[i][0]);
			d.sort_field(&m_pressureAccel[i][0]);
		}
	}
}


void SimulationDataPCISPHtest::emittedParticles(FluidModel* model, const unsigned int startIndex)
{
	// initialize kappa values for new particles
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	for (unsigned int j = startIndex; j < model->numActiveParticles(); j++)
	{
		m_predV[fluidModelIndex][j] = model->getVelocity(j);
		m_predX[fluidModelIndex][j] = model->getPosition(j);
		m_densityAdv[fluidModelIndex][j] = 0.0f;
		m_pressure[fluidModelIndex][j] = 0.0f;
		m_pressureAccel[fluidModelIndex][j].setZero();
	}
}