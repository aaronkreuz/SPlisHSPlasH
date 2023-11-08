#include "SimulationDataDFSPHvanilla.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"


using namespace SPH;

SimulationDataDFSPHvanilla::SimulationDataDFSPHvanilla(void) :
	m_factor(),
	m_pressure(),
	m_pressure_V(),
	m_source_term(),
	m_source_term_div(),
	m_pressureAccel(),
	m_density_adv(),
	m_aii()
{
}

SimulationDataDFSPHvanilla::~SimulationDataDFSPHvanilla(void)
{
	cleanup();
}


void SimulationDataDFSPHvanilla::init()
{
	Simulation* sim = Simulation::getCurrent();
	int nModels = Simulation::getCurrent()->numberOfFluidModels();
	int nBoundaries = Simulation::getCurrent()->numberOfBoundaryModels();

	m_factor.resize(nModels);
	m_density_adv.resize(nModels);
	m_pressure.resize(nModels);
	m_pressure_V.resize(nModels);
	m_source_term.resize(nModels);
	m_source_term_div.resize(nModels);
	m_pressureAccel.resize(nModels);
	m_aii.resize(nModels);

	for (auto i = 0; i < nModels; i++) {
		FluidModel* fm = Simulation::getCurrent()->getFluidModel(i);
		m_factor[i].resize(fm->numParticles(), 0.0);
		m_density_adv[i].resize(fm->numParticles(), 0.0);
		m_pressure[i].resize(fm->numParticles(), 0.0);
		m_pressure_V[i].resize(fm->numParticles(), 0.0);
		m_source_term[i].resize(fm->numParticles(), 0.0);
		m_source_term_div[i].resize(fm->numParticles(), 0.0);
		m_pressureAccel[i].resize(fm->numParticles(), Vector3r::Zero());
		m_aii[i].resize(fm->numParticles(), 0.0);
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
		m_pressure_V[i].clear();
		m_source_term[i].clear();
		m_source_term_div[i].clear();
		m_pressureAccel[i].clear();
		m_aii[i].clear();
	}

	m_factor.clear();
	m_density_adv.clear();
	m_pressure.clear();
	m_pressure_V.clear();
	m_source_term.clear();
	m_source_term_div.clear();
	m_pressureAccel.clear();
	m_aii.clear();
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
			m_pressure_V[i][j] = 0.0;
			m_factor[i][j] = 0.0;
			m_pressureAccel[i][j].setZero();
			m_source_term[i][j] = 0.0;
			m_source_term_div[i][j] = 0.0;
			m_aii[i][j] = 0.0;
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
			d.sort_field(&m_pressure_V[i][0]);
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
		m_pressure_V[fluidModelIndex][j] = 0.0;
	}
}
