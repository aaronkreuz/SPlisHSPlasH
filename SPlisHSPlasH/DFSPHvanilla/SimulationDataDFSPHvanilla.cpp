#include "SimulationDataDFSPHvanilla.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/Simulation.h"

using namespace SPH;

SimulationDataDFSPHvanilla::SimulationDataDFSPHvanilla(void) :
	m_factor(),
	m_pressure(),
	m_pressure_rho2_V(),
	m_pressureAccel(),
	m_density_adv()
{
}

SimulationDataDFSPHvanilla::~SimulationDataDFSPHvanilla(void)
{
	cleanup();
}


void SimulationDataDFSPHvanilla::init()
{
	int nModels = Simulation::getCurrent()->numberOfFluidModels();

	m_factor.resize(nModels);
	m_density_adv.resize(nModels);
	m_pressure.resize(nModels);
	m_pressure_rho2_V.resize(nModels);
	m_pressureAccel.resize(nModels);

	for (auto i = 0; i < nModels; i++) {
		FluidModel* fm = Simulation::getCurrent()->getFluidModel(i);
		m_factor[i].resize(fm->numParticles(), 0.0);
		m_density_adv[i].resize(fm->numParticles(), 0.0);
		m_pressure[i].resize(fm->numParticles(), 0.0);
		m_pressure_rho2_V[i].resize(fm->numParticles(), 0.0);
		m_pressureAccel[i].resize(fm->numParticles(), Vector3r::Zero());
	}
}

void SimulationDataDFSPHvanilla::cleanup()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		m_factor[i].clear();
		m_density_adv[i].clear();
		m_pressure[i].clear();
		m_pressure_rho2_V[i].clear();
		m_pressureAccel[i].clear();
	}
	m_factor.clear();
	m_density_adv.clear();
	m_pressure.clear();
	m_pressure_rho2_V.clear();
	m_pressureAccel.clear();
}

void SimulationDataDFSPHvanilla::reset()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel* fm = sim->getFluidModel(i);
		for (unsigned int j = 0; j < fm->numActiveParticles(); j++)
		{
			m_density_adv[i][j] = 0.0;
			m_pressure[i][j] = 0.0;
			m_pressure_rho2_V[i][j] = 0.0;
			m_factor[i][j] = 0.0;
			m_pressureAccel[i][j].setZero();
		}
	}
}

void SimulationDataDFSPHvanilla::performNeighborhoodSearchSort()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int i = 0; i < nModels; i++)
	{
		FluidModel* fm = sim->getFluidModel(i);
		const unsigned int numPart = fm->numActiveParticles();
		if (numPart != 0)
		{
			auto const& d = sim->getNeighborhoodSearch()->point_set(fm->getPointSetIndex());
			//d.sort_field(&m_factor[i][0]);
			//d.sort_field(&m_density_adv[i][0]);
			d.sort_field(&m_pressure[i][0]);
			d.sort_field(&m_pressure_rho2_V[i][0]);
		}
	}
}

void SimulationDataDFSPHvanilla::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	// initialize kappa values for new particles
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	for (unsigned int j = startIndex; j < model->numActiveParticles(); j++)
	{
		m_pressure[fluidModelIndex][j] = 0.0;
		m_pressure_rho2_V[fluidModelIndex][j] = 0.0;
	}
}
