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
		model->addField({ "p / rho^2", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureRho2(fluidModelIndex, i); }, true });
		model->addField({ "p_v / rho^2", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureRho2_V(fluidModelIndex, i); }, true });
		model->addField({ "pressure acceleration", FieldType::Vector3, [this, fluidModelIndex](const unsigned int i) -> Real* { return &m_simulationData.getPressureAccel(fluidModelIndex, i)[0]; } });
	}

	// Precompute values
	//////////////////////////////////////////////////////////////////////////
	TimeManager* tm = TimeManager::getCurrent();
	const Real h = tm->getTimeStepSize();

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

	performNeighborhoodSearch();

#ifdef USE_PERFORMANCE_OPTIMIZATION
	//////////////////////////////////////////////////////////////////////////
	// precompute the values V_j * grad W_ij for all neighbors
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("precomputeValues")
		precomputeValues();
	STOP_TIMING_AVG
#endif

	//////////////////////////////////////////////////////////////////////////
	// compute volume/density maps boundary contribution
	 if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	 	computeVolumeAndBoundaryX();
	 else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	 	computeDensityAndGradient();
	/////////////////////////////////////////////////////////////////////////

	// compute densities
	for (auto fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
		// density stored in fluid model
		computeDensities(fluidModelIndex);
	}

	// compute alpha factors required for stiffness parameter computation
	for (auto fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
		computeDFSPHFactor(fluidModelIndex);
	}

	// TODO: DIVERGENCE SOLVER

	// Reset accelerations and add gravity
	for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
		clearAccelerations(fluidModelIndex);

	// compute non-pressure forces
	sim->computeNonPressureForces();

	// update time step size
	sim->updateTimeStepSize();

	// advect velocities
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

	//////////////////////////////////////////////////////////////////////////
	// Constant Density Solver
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("constant density solver");
	pressureSolve();
	STOP_TIMING_AVG
	//////////////////////////////////////////////////////////////////////////
	// End Constant Density Solver
	//////////////////////////////////////////////////////////////////////////

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
}

void TimeStepDFSPHvanilla::pressureSolve()
{
	Simulation* sim = Simulation::getCurrent();
	TimeManager* tm = TimeManager::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const Real h = tm->getTimeStepSize();
	const Real invH = 1.0 / h;
	const Real invH2 = 1.0 / (h * h);


	// precompute advected densities
	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
		FluidModel* fm = sim->getFluidModel(fluidModelIndex);
		const unsigned int numParticles = fm->numActiveParticles();

		for (int i = 0; i < numParticles; i++) {
            computeDensityAdv(fluidModelIndex, i, h, fm->getDensity0());

			// no warm start
			m_simulationData.getPressureRho2(fluidModelIndex, i) = 0.0;
        }
	}

	m_iterations = 0;

	Real avg_density_err = 0.0;
	Real density_err = 0.0;
	bool chk = false;

	// loop
	while ((!chk || (m_iterations < m_minIterations)) && (m_iterations < m_maxIterations)) {

		chk = true;

		// Parallelize
		for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {
			FluidModel* fm = sim->getFluidModel(fluidModelIndex);
			const unsigned int numParticles = fm->numActiveParticles();

			// compute pressure accelerations
			for (int i = 0; i < numParticles; i++) {
                computePressureAccel(fluidModelIndex, i, fm->getDensity0(), m_simulationData.getPressureRho2Data());
            }
		}

		for (unsigned int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++)
		{
			avg_density_err = 0.0;

			FluidModel* fm = sim->getFluidModel(fluidModelIndex);
			Real density_err = 0.0;

			Real rho0 = fm->getDensity0();

			for (int i = 0; i < fm->numActiveParticles(); i++) {
				const Real& densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
				// determine source term
				const Real s_i = static_cast<Real>(1.0) - densityAdv;
				// compute matrix element
				Real aij_pj = compute_aij_pj(fluidModelIndex, i);
				aij_pj *= h * h;
				
				Real& p_rho2_i = m_simulationData.getPressureRho2(fluidModelIndex, i);
				const Real residuum = min(s_i - aij_pj, static_cast<Real>(0.0));     // r = b - A*p
				//p_rho2_i -= residuum * m_simulationData.getFactor(fluidModelIndex, i);

				p_rho2_i = max(p_rho2_i - 0.5 * (s_i - aij_pj) * m_simulationData.getFactor(fluidModelIndex, i), 0.0);

				//////////////////////////////////////////////////////////////////////////
				// Compute the sum of the density errors
				//////////////////////////////////////////////////////////////////////////
				density_err -= rho0 * residuum;
			}

			avg_density_err = density_err / fm->numActiveParticles();

			const Real eta = m_maxError * static_cast<Real>(0.01) * rho0;
			chk = chk && (avg_density_err <= eta);

		}

		m_iterations++;
	}

	// TODO: omp
	for (int fluidModelIndex = 0; fluidModelIndex < nFluids; fluidModelIndex++) {

		FluidModel* fm = sim->getFluidModel(fluidModelIndex);
        const unsigned int numParticles = fm->numActiveParticles();
		const Real rho0 = fm->getDensity0();

		for (int i = 0; i < numParticles; i++) {
			computePressureAccel(fluidModelIndex, i, rho0, m_simulationData.getPressureRho2Data(), true);
			fm->getVelocity(i) += h * m_simulationData.getPressureAccel(fluidModelIndex, i);
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

// void computeDFSPHFactor(const unsigned int fluidModelIndex)
// {
// 	Simulation* sim = Simulation::getCurrent();
// 	std::vector<Real>& factors_model = m_simulationData.getFactorsModel(fluidModelIndex);
// 	FluidModel* model = sim->getFluidModel(fluidModelIndex);
//     const unsigned int numParticles = model->numActiveParticles();
// 	const unsigned int nFluids = sim->numberOfFluidModels();
// 
// 	for (unsigned int i = 0; i < numParticles; i++) {
// 		Real& alpha_i = factors_model[i];
// 		const Vector3r xi = model->getPosition(i);
// 
// 		const Real density_i = model->getDensity(i);
// 
// 		Vector3r sumGradKernel = Vector3r::Zero();
// 		Real sumGradKernel2 = static_cast<Real>(0.0);
// 
// 		forall_fluid_neighbors(
// 			const Vector3r massGradKernel = model->getMass(neighborIndex) * sim->gradW(xi - xj);
// 			sumGradKernel += massGradKernel;
// 			sumGradKernel2 += massGradKernel.squaredNorm();
// 		);
// 
// 		// TODO: Boundary neighbors
// 
// 		Real sumGradKernelSqrNorm = sumGradKernel.squaredNorm();
// 
// 		// TODO: Maybe add small constant to avoid numerical issues to denominator
// 		alpha_i = density_i / (sumGradKernelSqrNorm + sumGradKernel2);
// 	}
// 
// }

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

			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					grad_p_i -= gradRho;
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					const Vector3r grad_p_j = -Vj * sim->gradW(xi - xj);
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

/** Compute rho_adv,i^(0) (see equation in Section 3.3 in [BK17])
  * using the velocities after the non-pressure forces were applied.
**/
void TimeStepDFSPHvanilla::computeDensityAdv(const unsigned int fluidModelIndex, const unsigned int i, const Real h, const Real density0)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const Real& density = model->getDensity(i);
	Real& densityAdv = m_simulationData.getDensityAdv(fluidModelIndex, i);
	const Vector3r& xi = model->getPosition(i);
	const Vector3r& vi = model->getVelocity(i);
	Real delta = 0.0;
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors(
		const Vector3r & vj = fm_neighbor->getVelocity(neighborIndex);
	delta += (vi - vj).dot(sim->gradW(xi - xj));
	//delta += fm_neighbor->getVolume(neighborIndex) * (vi - vj).dot(sim->gradW(xi - xj));
	);
	// assumes that all fluid particles have the same volume
	delta *= model->getVolume(i);

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors(
			const Vector3r & vj = bm_neighbor->getVelocity(neighborIndex);
		delta += bm_neighbor->getVolume(neighborIndex) * (vi - vj).dot(sim->gradW(xi - xj));
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	{
		forall_density_maps(
			Vector3r vj;
		bm_neighbor->getPointVelocity(xi, vj);
		delta -= (vi - vj).dot(gradRho);
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	{
		forall_volume_maps(
			Vector3r vj;
		bm_neighbor->getPointVelocity(xj, vj);
		delta += Vj * (vi - vj).dot(sim->gradW(xi - xj));
		);
	}

	densityAdv = density / density0 + h * delta;
}

// void TimeStepDFSPHvanilla::computeDensityAdv(const unsigned int fluidModelIndex, const unsigned int i, const Real h, const Real density0)
// {
// 	Simulation* sim = Simulation::getCurrent();
// 	const unsigned int nFluids = sim->numberOfFluidModels();
// 	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
// 
// 	FluidModel* model = Simulation::getCurrent()->getFluidModel(fluidModelIndex);
// 	Real& density_adv = m_simulationData.getDensityAdv(fluidModelIndex, i);
// 	const Real density = model->getDensity(i);
// 
// 	const Vector3r xi = model->getPosition(i);
// 
// 	Real deltaDensity = 0.0;
// 
// 	// not avx optimized
// 	forall_fluid_neighbors(
// 		Real mass_j = model->getMass(neighborIndex);
// 		const Vector3r diffPredV = model->getVelocity(neighborIndex) - model->getVelocity(i);
// 		deltaDensity += mass_j * diffPredV.dot(sim->gradW(xi - xj));
// 	);
// 
// 	//////////////////////////////////////////////////////////////////////////
// 	// Boundary
// 	 if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
// 	 {
// 	 	forall_boundary_neighbors(
// 			const Real V_j = bm_neighbor->getVolume(neighborIndex);
// 	 		const Vector3r v_j = bm_neighbor->getVelocity(neighborIndex);
// 	 		deltaDensity += V_j * (model->getVelocity(i) - v_j).dot(CubicKernel::gradW(xi - xj));
// 	 	);
// 	 }
// 	 
// 	 if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
// 	 {
// 	 	forall_density_maps(
// 	 		Vector3r vj;
// 	 		bm_neighbor->getPointVelocity(xi, vj);
// 	 		deltaDensity -= (model->getVelocity(i) - vj).dot(gradRho);
// 	 	);
// 	 }
// 	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
// 	{
// 		forall_volume_maps(
// 			Vector3r vj;
// 			bm_neighbor->getPointVelocity(xj, vj);
// 			deltaDensity += Vj * (model->getVelocity(i) - vj).dot(sim->gradW(xi - xj));
// 		);
// 	}
// 	//////////////////////////////////////////////////////////////////////////
// 
// 	density_adv = density + h * deltaDensity;
// }

void TimeStepDFSPHvanilla::computeDensityChange(const unsigned int fluidModelIndex, const unsigned int i, const Real h)
{
}


/** Compute pressure accelerations using the current pressure values of the particles
 */
void TimeStepDFSPHvanilla::computePressureAccel(const unsigned int fluidModelIndex, const unsigned int i, const Real density0, std::vector<std::vector<Real>>& pressure_rho2, const bool applyBoundaryForces)
{
	//////////////////////////////////////////////////////////////////////////
	// Sum_j [ m_j*(kappa_i + kappa_j)*gradW_ij ]
	//////////////////////////////////////////////////////////////////////////

	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
    const unsigned int numParticles = model->numActiveParticles();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	Vector3r& ap_i = m_simulationData.getPressureAccel(fluidModelIndex, i);
	ap_i.setZero();

	const Real& kappa_i = m_simulationData.getPressureRho2(fluidModelIndex, i);
	const Vector3r& xi = model->getPosition(i);

	Vector3r sum = Vector3r::Zero();
	
	// parallelize!
	forall_fluid_neighbors(
		const Real& m_j = model->getMass(neighborIndex);
		const Real& kappa_j = m_simulationData.getPressureRho2(pid, neighborIndex);
		const Vector3r gradW = sim->gradW(xi - xj);

		ap_i -= m_j * (kappa_i + kappa_j) * gradW;
	);

	//////////////////////////////////////////////////////////////////////////
	// Boundary Handling
	
	if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
	{
		forall_boundary_neighbors(
			const Vector3r grad_p_j = -bm_neighbor->getVolume(neighborIndex) * sim->gradW(xi - xj);

			const Vector3r a = (Real)1.0 * kappa_i * grad_p_j;
			ap_i += a;
			if (applyBoundaryForces)
				bm_neighbor->addForce(xj, -model->getMass(i) * a);
		);
	}	
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	{
		forall_density_maps(
			const Vector3r a = (Real)1.0 * kappa_i * gradRho;
			ap_i += a;
			if (applyBoundaryForces)
				bm_neighbor->addForce(xj, -model->getMass(i) * a);
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	{
		forall_volume_maps(
			const Vector3r grad_p_j = -Vj * sim->gradW(xi - xj);
			const Vector3r a = (Real)1.0 * kappa_i * grad_p_j;
			ap_i += a;

			if (applyBoundaryForces)
				bm_neighbor->addForce(xj, -model->getMass(i) * a);
		);
	}
    
	
	//////////////////////////////////////////////////////////////////////////
}

// original
Real TimeStepDFSPHvanilla::compute_aij_pj(const unsigned int fluidModelIndex, const unsigned int i)
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = sim->getFluidModel(fluidModelIndex);
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

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
		const Vector3r & aj = m_simulationData.getPressureAccel(pid, neighborIndex);
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
			aij_pj += bm_neighbor->getVolume(neighborIndex) * ai.dot(sim->gradW(xi - xj));
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
	{
		forall_density_maps(
			aij_pj -= ai.dot(gradRho);
		);
	}
	else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
	{
		forall_volume_maps(
			aij_pj += Vj * ai.dot(sim->gradW(xi - xj));
		);
	}
	return aij_pj;
}

