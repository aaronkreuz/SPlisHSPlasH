#include "TimeStepDFSPHvanilla.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/SPHKernels.h"
#include "SimulationDataDFSPHvanilla.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"


using namespace SPH;
using namespace std;
using namespace GenParam;


int TimeStepDFSPHvanilla::SOLVER_ITERATIONS_V = -1;
int TimeStepDFSPHvanilla::MAX_ITERATIONS_V = -1;
int TimeStepDFSPHvanilla::MAX_ERROR_V = -1;
int TimeStepDFSPHvanilla::USE_DIVERGENCE_SOLVER = -1;


TimeStepDFSPHvanilla::TimeStepDFSPHvanilla() :
	TimeStep(),
	m_simulationData()
{
	m_simulationData.init();

	m_counter = 0;
	m_iterationsV = 0;
	m_enableDivergenceSolver = true;
	m_maxIterationsV = 100;
	m_maxErrorV = static_cast<Real>(0.1);

	Simulation * sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		model->addField({ "factor", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getFactor(fluidModelIndex, i); } });
		model->addField({ "advected density", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getDensityAdv(fluidModelIndex, i); } });
		model->addField({ "p / rho^2", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressure(fluidModelIndex, i); }, true });
		model->addField({ "p_v / rho^2", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureRho2_V(fluidModelIndex, i); }, true });
		model->addField({ "pressure acceleration", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureAccel(fluidModelIndex, i)[0]; } });
	}

	// Precompute values
	//////////////////////////////////////////////////////////////////////////
	TimeManager* tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	performNeighborhoodSearch();

	for (int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++){
		computeDensities(fluidModelIndex);
	}

	// compute diagonal matrix elements for t0
	for (int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++){
		for (int i = 0; i < Simulation::getCurrent()->getFluidModel(fluidModelIndex)->numActiveParticles(); i++) {
			compute_aii(fluidModelIndex, i, h);
		}
	}
	//////////////////////////////////////////////////////////////////////////
}

TimeStepDFSPHvanilla::~TimeStepDFSPHvanilla(void)
{
	// remove all particle fields
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++)
	{
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		model->removeFieldByName("factor");
		model->removeFieldByName("advected density");
		model->removeFieldByName("p / rho^2");
		model->removeFieldByName("p_v / rho^2");
		model->removeFieldByName("pressure acceleration");
	}
}

void TimeStepDFSPHvanilla::initParameters()
{
	// TODO: add parameters
}

void TimeStepDFSPHvanilla::step()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	TimeManager* tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	// reset the accelerations of the particles and add gravity
	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
		clearAccelerations(fluidModelIndex);
	}

	// compute non-pressure forces
	// INFO: stored in acceleration array of fluid-model
	// note that neighbors and densities are already determined at this point
	// only viscosity (XSPH)
	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
		for (int i = 0; i < Simulation::getCurrent()->getFluidModel(fluidModelIndex)->numActiveParticles(); i++){
			computeViscosityForce(fluidModelIndex, i, h);
		}
	}

	// advect velocities based on non-pressure forces
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = model->numActiveParticles();

		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)
			for (int i = 0; i < numParticles; i++)
			{
				if (model->getParticleState(i) == ParticleState::Active)
				{
					Vector3r& vel = model->getVelocity(i);
					// accelerations already containing gravity and non-pressure accels
					vel += h * model->getAcceleration(i);
				}
			}
		}
	}

	// compute constant density source term [(rho0 - rho⁺) / h]

	// compute advected densities after applying non-pressure forces [rho^*]
	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
		FluidModel* fm = Simulation::getCurrent()->getFluidModel(fluidModelIndex);
		for (int i = 0; i < fm->numActiveParticles(); i++) {
			computeDensityAdv(fluidModelIndex, i, h, fm->getDensity0());

			// no warm start
			m_simulationData.getPressure(fluidModelIndex, i) = 0.0;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Constant Density Solver
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("constant density solver");
	pressureSolve();
	STOP_TIMING_AVG;


	// update velocities using pressure accelerations
	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
		FluidModel* fm = Simulation::getCurrent()->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = fm->numActiveParticles();

		for (int i = 0; i < numParticles; i++) {
			if (fm->getParticleState(i) == ParticleState::Active)
			{
				fm->getVelocity(i) += h * m_simulationData.getPressureAccel(fluidModelIndex, i);
			}
		}
	}

	// update positions
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
		FluidModel* model = sim->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = model->numActiveParticles();

		for (unsigned int i = 0; i < numParticles; i++) {
			if (model->getParticleState(i) == ParticleState::Active)
			{
				Vector3r &xi = model->getPosition(i);
				const Vector3r &vi = model->getVelocity(i);
				xi += h * vi;
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// emit new particles and perform an animation field step
	//////////////////////////////////////////////////////////////////////////
	sim->emitParticles();
	sim->animateParticles();

	//////////////////////////////////////////////////////////////////////////
	// Compute new time
	//////////////////////////////////////////////////////////////////////////
	tm->setTime(tm->getTime() + h);

	performNeighborhoodSearch();

	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
		computeDensities(fluidModelIndex);
	}

	for(int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
		FluidModel* fm = sim->getFluidModel(fluidModelIndex);
		for (int i = 0; i < fm->numActiveParticles(); i++){
			compute_aii(fluidModelIndex, i, h);
		}
	}

	// TODO: DIVERGENCE SOLVER

	// //////////////////////////////////////////////////////////////////////////
	// // compute volume/density maps boundary contribution
	//  if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	//  	computeVolumeAndBoundaryX();
	//  else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	//  	computeDensityAndGradient();
	// /////////////////////////////////////////////////////////////////////////

	// // compute densities
	// for (auto fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
	// 	// density stored in fluid model
	// 	computeDensities(fluidModelIndex);
	// }

	// // compute alpha factors required for stiffness parameter computation
	// for (auto fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
	// 	computeDFSPHFactor(fluidModelIndex);
	// }


	// // update time step size
	// sim->updateTimeStepSize();

}

void TimeStepDFSPHvanilla::computeViscosityForce(const unsigned int fluidModelIndex, const unsigned i, const Real h) {
	Simulation* sim = Simulation::getCurrent();
	FluidModel* fm = Simulation::getCurrent()->getFluidModel(fluidModelIndex);
	Vector3r& accel_i = fm->getAcceleration(i);
	auto nFluids = Simulation::getCurrent()->numberOfFluidModels();
	const Vector3r &xi = fm->getPosition(i);

	Real factor = 0.01 / h;
	Vector3r sum = Vector3r::Zero();

	forall_fluid_neighbors(
		sum += (fm_neighbor->getMass(neighborIndex) / fm_neighbor->getDensity(neighborIndex)) * (fm_neighbor->getVelocity(neighborIndex) - fm->getVelocity(i)) * CubicKernel::W(xi - xj);
	);

	accel_i += factor * sum;
}

void TimeStepDFSPHvanilla::pressureSolve(){
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	TimeManager* tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

	bool solved = false;
	m_iterations = 0;

	Real density_err = 0.0;
	Real avg_density_err = 0.0;

	while ((!solved || (m_iterations < m_minIterations)) && m_iterations < m_maxIterations){
		solved = true;
		avg_density_err = 0.0;

		// compute pressure accelerations
		for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
			FluidModel* fm = sim->getFluidModel(fluidModelIndex);
			const unsigned int numParticles = fm->numActiveParticles();

			for (int i = 0; i < numParticles; i++) {
				computePressureAccel(fluidModelIndex, i, fm->getDensity0(), m_simulationData.getPressureRho2Data());
			}
		}

		 for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
			FluidModel* fm = sim->getFluidModel(fluidModelIndex);
			const unsigned int numParticles = fm->numActiveParticles();
			Real density0 = fm->getDensity0();

			density_err = 0.0;
			for(int i = 0; i < numParticles; i++) {
				Real density_adv = m_simulationData.getDensityAdv(fluidModelIndex, i);
				Real& pressure_i = m_simulationData.getPressure(fluidModelIndex, i);
				Real& diag_i = m_simulationData.getDiagElement(fluidModelIndex, i);
				Real aij_pj = compute_aij_pj(fluidModelIndex, i);

				if (abs(m_simulationData.getDiagElement(fluidModelIndex, i)) > 1.0e-6){
					pressure_i += (0.5 / diag_i) * (/*source term*/ ((density0 - density_adv) / h) - aij_pj);
				}
				else {
					pressure_i = 0.0;
				}

				// clamp to zero
				pressure_i = max(pressure_i, 0.0f);
				density_err += max(/*source term*/ ((density0 - density_adv) / h) - aij_pj, 0.0f) * h;
			}

			avg_density_err = density_err / numParticles;
			solved = solved && (avg_density_err <= (m_maxError * 0.01 * density0));
		 }

		m_iterations++;
	}

	// final acceleration update
	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++){
		FluidModel* fm = sim->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = fm->numActiveParticles();

		for (int i = 0; i < numParticles; i++) {
			computePressureAccel(fluidModelIndex, i, fm->getDensity0(), m_simulationData.getPressureRho2Data());
		}
	}
}

void TimeStepDFSPHvanilla::pressureSolveIteration(const unsigned int fluidModelIndex, Real& avg_density_err)
{
}

void TimeStepDFSPHvanilla::divergenceSolve()
{
}

void TimeStepDFSPHvanilla::divergenceSolveIteration(const unsigned int fluidModelIndex, Real &avg_density_err)
{
}

void TimeStepDFSPHvanilla::reset()
{
	TimeStep::reset();
	m_simulationData.reset();
	m_counter = 0;
	m_iterations = 0;
	m_iterationsV = 0;
}

void TimeStepDFSPHvanilla::performNeighborhoodSearch()
{
	if (Simulation::getCurrent()->zSortEnabled())
	{
		if (m_counter % 500 == 0)
		{
			Simulation::getCurrent()->performNeighborhoodSearchSort();
			m_simulationData.performNeighborhoodSearchSort();
		}
		m_counter++;
	}

	Simulation::getCurrent()->performNeighborhoodSearch();
}

void TimeStepDFSPHvanilla::emittedParticles(FluidModel *model, const unsigned int startIndex)
{
	m_simulationData.emittedParticles(model, startIndex);
}

void TimeStepDFSPHvanilla::resize()
{
	m_simulationData.init();
}

void TimeStepDFSPHvanilla::computeDFSPHFactor(const unsigned int fluidModelIndex)
{
	//////////////////////////////////////////////////////////////////////////
	// Init parameters
	//////////////////////////////////////////////////////////////////////////

	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const int numParticles = (int)model->numActiveParticles();

#pragma omp parallel default(shared)
	{
		//////////////////////////////////////////////////////////////////////////
		// Compute pressure stiffness denominator
		//////////////////////////////////////////////////////////////////////////

#pragma omp for schedule(static)  
		for (int i = 0; i < numParticles; i++)
		{
			//////////////////////////////////////////////////////////////////////////
			// Compute gradient dp_i/dx_j * (1/kappa)  and dp_j/dx_j * (1/kappa)
			// (see Equation (8) and the previous one [BK17])
			// Note: That in all quantities rho0 is missing due to our
			// implementation of multiphase simulations.
			//////////////////////////////////////////////////////////////////////////
			const Vector3r& xi = model->getPosition(i);
			Real sum_grad_p_k = 0.0;
			Vector3r grad_p_i;
			grad_p_i.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				const Vector3r grad_p_j = -fm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
			sum_grad_p_k += grad_p_j.squaredNorm();
			grad_p_i -= grad_p_j;
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					const Vector3r grad_p_j = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);
				grad_p_i -= grad_p_j;
				);
			}

			sum_grad_p_k += grad_p_i.squaredNorm();

			//////////////////////////////////////////////////////////////////////////
			// Compute factor as: factor_i = -1 / (a_ii * rho_i^2)
			// where a_ii is the diagonal entry of the linear system 
			// for the pressure A * p = source term
			//////////////////////////////////////////////////////////////////////////
			Real& factor = m_simulationData.getFactor(fluidModelIndex, i);
			if (sum_grad_p_k > m_eps)
				factor = static_cast<Real>(1.0) / (sum_grad_p_k);
			else
				factor = 0.0;
		}
	}
}


void TimeStepDFSPHvanilla::computeDensityAdv(const unsigned int fluidModelIndex, const unsigned int i, const Real h, const Real density0)
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	FluidModel* model = Simulation::getCurrent()->getFluidModel(fluidModelIndex);
	const Real rho0 = model->getDensity0();
	Real& density_adv = m_simulationData.getDensityAdv(fluidModelIndex, i);
	const Real density = model->getDensity(i);

	const Vector3r xi = model->getPosition(i);

	Real deltaDensity = 0.0;

	forall_fluid_neighbors(
		// Real V_j = model->getVolume(neighborIndex);
		Real mass = fm_neighbor->getMass(neighborIndex);
		const Vector3r diffPredV = fm_neighbor->getVelocity(neighborIndex) - model->getVelocity(i);
		deltaDensity += mass * diffPredV.dot(sim->gradW(xi - xj));
	);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	 if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	 {
		forall_boundary_neighbors(
			const Real V_j = bm_neighbor->getVolume(neighborIndex);
			const Real mass_j = V_j * rho0;
			const Vector3r v_j = bm_neighbor->getVelocity(neighborIndex);
			deltaDensity += mass_j * (model->getVelocity(i) - v_j).dot(CubicKernel::gradW(xi - xj));
		);
	 }
	//////////////////////////////////////////////////////////////////////////

	density_adv = density + h * deltaDensity;
}

void TimeStepDFSPHvanilla::computeDensityChange(const unsigned int fluidModelIndex, const unsigned int i, const Real h)
{
}


/** Compute pressure accelerations using the current pressure values of the particles
 */
void TimeStepDFSPHvanilla::computePressureAccel(const unsigned int fluidModelIndex, const unsigned int i, const Real density0, std::vector<std::vector<Real>>& pressure_rho2, const bool applyBoundaryForces)
{
	//////////////////////////////////////////////////////////////////////////
	// Sum_j [ m_j*(p_i/rho0_i² + p_j/rho0_j²)*gradW_ij ]
	//////////////////////////////////////////////////////////////////////////

	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int nFluids = sim->numberOfFluidModels();

	Vector3r& ap_i = m_simulationData.getPressureAccel(fluidModelIndex, i);
	ap_i.setZero();

	const Real& p_i = m_simulationData.getPressure(fluidModelIndex, i);
	const Real dpi = p_i / (model->getDensity(i) * model->getDensity(i));

	const Vector3r& xi = model->getPosition(i);

	Vector3r sum = Vector3r::Zero();
	
	forall_fluid_neighbors(
		const Real& V_j = model->getVolume(neighborIndex);
		const Real& mass = model->getMass(neighborIndex);
		const Real& p_j = m_simulationData.getPressure(pid, neighborIndex);
		const Real dpj = p_j / (fm_neighbor->getDensity(neighborIndex) * fm_neighbor->getDensity(neighborIndex));
		const Vector3r gradW = sim->gradW(xi - xj);

		ap_i -= mass * (dpi + dpj) * gradW;
	);

	//////////////////////////////////////////////////////////////////////////
	// Boundary Handling
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors(
			const Vector3r grad_p_j = bm_neighbor->getVolume(neighborIndex) * density0 * sim->gradW(xi - xj);

			const Vector3r a = (Real)1.0 * dpi * grad_p_j;
			ap_i -= a;
			if (applyBoundaryForces)
				bm_neighbor->addForce(xj, -model->getMass(i) * a);
		);
	}	
	//////////////////////////////////////////////////////////////////////////
	// ap_i *= model->getDensity0();
}

void TimeStepDFSPHvanilla::compute_aii(const unsigned int fluidModelIndex, const unsigned int i, const Real h){
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int nFluids = sim->numberOfFluidModels();

	const Real& rho_i = model->getDensity(i);
	const Real& rho0_i = model->getDensity0();
	const Vector3r& xi = model->getPosition(i);

	Real factor = -(h *rho0_i ) / (rho_i * rho_i);
	Real sum1 = 0.0;
	Vector3r sum2 = Vector3r::Zero();

	forall_fluid_neighbors(
		sum1 += (fm_neighbor->getMass(neighborIndex) * sim->gradW(xi-xj)).squaredNorm();
		sum2 += fm_neighbor->getMass(neighborIndex) * sim->gradW(xi-xj);
	);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors(
			const Vector3r grad_p_j = bm_neighbor->getVolume(neighborIndex)*rho0_i * sim->gradW(xi - xj);
			sum2 += grad_p_j;
		);
	}
	//////////////////////////////////////////////////////////////////////////

	m_simulationData.setDiagElement(fluidModelIndex, i, factor * (sum1 + sum2.squaredNorm()));
}

// original
Real TimeStepDFSPHvanilla::compute_aij_pj(const unsigned int fluidModelIndex, const unsigned int i)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const Real density0 = model->getDensity0();

	//////////////////////////////////////////////////////////////////////////
	// Compute A*p which is the change of the density when applying the 
	// pressure forces. 
	// \sum_j a_ij * p_j = h^2 \sum_j V_j (a_i - a_j) * gradW_ij
	// This is the RHS of Equation (12) in [BK17]
	//////////////////////////////////////////////////////////////////////////
	const Vector3r& xi = model->getPosition(i);
	const Vector3r& ai = m_simulationData.getPressureAccel(fluidModelIndex, i);
	Real aij_pj = 0.0;

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors(
		const Vector3r& aj = m_simulationData.getPressureAccel(pid, neighborIndex);
		//aij_pj += fm_neighbor->getVolume(neighborIndex) * (ai - aj).dot(sim->gradW(xi - xj));
		aij_pj += (ai - aj).dot(sim->gradW(xi - xj));
	);
	// assumes that all fluid particles have the same volume
	aij_pj *= model->getVolume(i);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors(
			aij_pj += bm_neighbor->getVolume(neighborIndex)* density0 * ai.dot(sim->gradW(xi - xj));
		);
	}

	return aij_pj;
}

