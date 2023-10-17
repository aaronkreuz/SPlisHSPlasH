#include "TimeStepPCISPHtest.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataPCISPHtest.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "../Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"

using namespace SPH;
using namespace std;
using namespace GenParam;

int TimeStepPCISPHtest::ENABLE_DEBUG_OUTPUT = -1;

TimeStepPCISPHtest::TimeStepPCISPHtest() :
	TimeStep()
{
	m_simulationData.init();
	std::cout << "init simulation data done" << std::endl;
	m_counter = 0;
	m_minIterations = 3;

	Simulation *sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		
		// INFO: addField (FieldDescription::name, FieldDescription::type, FieldDescription::getFunction)
		// INFO: FieldDescription::getFunction returns a pointer to the function that returns the value of the field
	
		model->addField({ "density_adv", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getDensityAdv(fluidModelIndex, i); } });
		model->addField({ "pressure", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressure(fluidModelIndex, i); } });
		model->addField({ "pressure acceleration", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureAccel(fluidModelIndex, i)[0]; } });
	}
}

TimeStepPCISPHtest::~TimeStepPCISPHtest(void)
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);

		model->removeFieldByName("density_adv");
		model->removeFieldByName("pressure");
		model->removeFieldByName("pressure acceleration");
	}
}

void TimeStepPCISPHtest::step()
{
	// Get required classes
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	TimeManager* tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	// Q: how is the actual neighborhood search performed?
	// 1) compute current neighborhood for particles
	performNeighborhoodSearch();

	// 2) reconstruct densities at particle positions for every fluid model
	for (auto i = 0; i < nFluids; i++) {
		// clear accelerations already adds gravity
		clearAccelerations(i);
		// densitites are stored in corresponding fluid model
		computeDensities(i);
	}

	// 3) compute non-pressure forces
	// stored in the acceleration field of the corresponding fluid model
	sim->computeNonPressureForces();

	// TODO: 3.1) already update velocity with nonPressureForce accelerations

	// TODO: 3.2) set pressure and pressureAcceleration to 0
	for(auto i = 0; i < nFluids; i++) {
		FluidModel* fModel = sim->getFluidModel(i);
		for(auto j = 0; j < (int)fModel->numParticles(); j++) {
			m_simulationData.getPressure(i, j) = 0.0;
			m_simulationData.getPressureAccel(i, j).setZero();
		}
	}

	m_iterations = 0;
	// 4) pressure solver
	bool pressureSolved = false;
	while (m_iterations < m_maxIterations && ((m_iterations < m_minIterations) || !pressureSolved)) {

		// 4.1) predict velocities and positions
		for (auto fluidIndex = 0; fluidIndex < nFluids; fluidIndex++) {
			FluidModel* fModel = sim->getFluidModel(fluidIndex);
			for (auto j = 0; j < (int)fModel->numParticles(); j++) {
				auto& predV = m_simulationData.getPredictedVelocity(fluidIndex, j);
				auto& predX = m_simulationData.getPredictedPosition(fluidIndex, j);
				const auto& velocity = fModel->getVelocity(j);
				const auto& position = fModel->getPosition(j);

				predV = velocity + h * (fModel->getAcceleration(j) + m_simulationData.getPressureAccel(fluidIndex, j));
				predX = position + h * predV;
            }
        }

		// TODO: Recompute time step size here? CFL or other method?

		// 4.2) Recompute smoothing kernel values


		Real error_sum = 0.0;

        // 4.3) pressure solver iteration
		for (auto fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
			FluidModel* fm = sim->getFluidModel(fluidModelIndex);
			for (auto i = 0; i < fm->numParticles(); i++) {
				// 4.3a) Compute predicted density
				// INFO: Predicted (advected) density is stored in the simulation data
				auto& advDensity = m_simulationData.getDensityAdv(fluidModelIndex, i);

				// NOTE: Naming of surrounding indices is fixed. I.e. fluidModelIndex, i, nFluids...
				// INFO: pid is the index of the fluid model of neighor particle, neighborIndex is the index of the neighbor particle
				forall_fluid_neighbors(
					advDensity += fm_neighbor->getMass(neighborIndex) * sim->W(m_simulationData.getPredictedPosition(fluidModelIndex, i) - m_simulationData.getPredictedPosition(pid, neighborIndex));
				);

				// TODO: Boundary handling

				// add error to error sum (rho0 - advDensity)
				error_sum += fm->getDensity0() - advDensity;

				// 4.3b) Compute pressure
				// INFO: Pressure is stored in the simulation data
				auto& pressure = m_simulationData.getPressure(fluidModelIndex, i);
				pressure = m_simulationData.getPcisphFactor(fluidModelIndex, i) * (fm->getDensity0() - advDensity);

				// 4.3c) Compute pressure acceleration
				// INFO: Pressure acceleration is stored in the simulation data
				auto& pressureAccel = m_simulationData.getPressureAccel(fluidModelIndex, i);
				pressureAccel.setZero();
				const Real pressureAccelFactor = (fm->getMass(i) * static_cast<Real>(2.0) * pressure) / (fm->getDensity0() * fm->getDensity0());
				Vector3r sumGradW = { 0.0, 0.0, 0.0 };

				forall_fluid_neighbors(
					sumGradW += sim->gradW(m_simulationData.getPredictedPosition(pid, neighborIndex));
				);

				pressureAccel = -pressureAccelFactor * sumGradW;
			}
			// 4.3c) Average error (per fluid model)
			// INFO: error_sum is the sum of all errors of all particles of all fluid models
			// INFO: error_sum is divided by the number of particles of all fluid models
			Real avg_error = error_sum / fm->numParticles();
			pressureSolved = pressureSolved && false;

		}


	}

	// 5) Update velocities and positions
	for (auto fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
		FluidModel* fm = sim->getFluidModel(fluidModelIndex);
		for (auto i = 0; i < fm->numParticles(); i++) {

			auto& position = fm->getPosition(i);
			auto& velocity = fm->getVelocity(i);

			// TODO: Add again the non-pressure acceleration?
			// fm->getAcceleration() contains the non-pressure accelerations
			velocity += h * (m_simulationData.getPressureAccel(fluidModelIndex, i) + fm->getAcceleration(i));
			position += h * velocity;
		}
	}
	
}

/*
void timeStepWCSPHstep()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	TimeManager* tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	performNeighborhoodSearch();

#ifdef USE_PERFORMANCE_OPTIMIZATION
	precomputeValues();
#endif



	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
		computeVolumeAndBoundaryX();
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
		computeDensityAndGradient();

	// Compute accelerations: a(t)
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		clearAccelerations(fluidModelIndex);
		computeDensities(fluidModelIndex);
	}
	sim->computeNonPressureForces();


	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		const Real density0 = model->getDensity0();
#pragma omp parallel default(shared)
		{
#pragma omp for schedule(static)  
			for (int i = 0; i < (int)model->numActiveParticles(); i++)
			{
				Real& density = model->getDensity(i);
				density = max(density, density0);
				m_simulationData.getPressure(fluidModelIndex, i) = m_stiffness * (pow(density / density0, m_exponent) - static_cast<Real>(1.0));
			}
		}

		computePressureAccels(fluidModelIndex);
	}

	sim->updateTimeStepSize();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
#pragma omp parallel default(shared)
		{
#pragma omp for schedule(static) 
			for (int i = 0; i < (int)model->numActiveParticles(); i++)
			{
				if (model->getParticleState(i) == ParticleState::Active)
				{
					Vector3r& pos = model->getPosition(i);
					Vector3r& vel = model->getVelocity(i);
					Vector3r& accel = model->getAcceleration(i);
					accel += m_simulationData.getPressureAccel(fluidModelIndex, i);
					vel += accel * h;
					pos += vel * h;
				}
			}
		}
	}

	sim->emitParticles();
	sim->animateParticles();

	// Compute new time	
	tm->setTime(tm->getTime() + h);
}
*/


void TimeStepPCISPHtest::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
}

/*
void TimeStepPCISPHtest::computePressureAccels(const unsigned int fluidModelIndex)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const Real density0 = model->getDensity0();
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	// Compute pressure forces
#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r& xi = model->getPosition(i);
			const Real density_i = model->getDensity(i);

			Vector3r& ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
			ai.setZero();

			const Real dpi = m_simulationData.getPressure(fluidModelIndex, i) / (density_i * density_i);
			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				// Pressure 
				const Real density_j = fm_neighbor->getDensity(neighborIndex) * density0 / fm_neighbor->getDensity0();
			const Real dpj = m_simulationData.getPressure(pid, neighborIndex) / (density_j * density_j);
			ai -= density0 * fm_neighbor->getVolume(neighborIndex) * (dpi + dpj) * sim->gradW(xi - xj);
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			const Real dpj = m_simulationData.getPressure(fluidModelIndex, i) / (density0 * density0);
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					const Vector3r a = density0 * bm_neighbor->getVolume(neighborIndex) * (dpi + dpj) * sim->gradW(xi - xj);
				ai -= a;
				bm_neighbor->addForce(xj, model->getMass(i) * a);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					const Vector3r a = -density0 * (dpi + dpj) * gradRho;
				ai -= a;
				bm_neighbor->addForce(xj, model->getMass(i) * a);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					const Vector3r a = density0 * Vj * (dpi + dpj) * sim->gradW(xi - xj);
				ai -= a;
				bm_neighbor->addForce(xj, model->getMass(i) * a);
				);
			}
		}
	}
}
*/

void TimeStepPCISPHtest::performNeighborhoodSearch()
{
	if (Simulation::getCurrent()->zSortEnabled())
	{
		// sort every 1000 steps
		if (m_counter % 1000 == 0)
		{
			Simulation::getCurrent()->performNeighborhoodSearchSort();
			m_simulationData.performNeighborhoodSearchSort();
		}
		m_counter++;
	}

	Simulation::getCurrent()->performNeighborhoodSearch();
}

void TimeStepPCISPHtest::emittedParticles(FluidModel* model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepPCISPHtest::resize()
{
	m_simulationData.init();
}
