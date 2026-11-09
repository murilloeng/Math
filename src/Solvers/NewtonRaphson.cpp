//std
#include <cstdio>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Solvers/NewtonRaphson.hpp"

//x: state
//r: residue
//p: continuation parameter

//target system: r(x, p) = 0

//tangent on x: K(x, p) = -dr/dx(x, p)
//tangent on p: g(x, p) = +dr/dp(x, p)

namespace math
{
	namespace solvers
	{
		//solve
		void NewtonRaphson::compute_step(void)
		{
			m_dps = m_dp;
			m_status = false;
			memcpy(m_dxs, m_dx, m_size * sizeof(double));
			for(m_attempt = 0; m_attempt < m_attempt_max; m_attempt++)
			{
				if(!predictor()) break;
				for(m_iteration = 0; m_iteration < m_iteration_max; m_iteration++)
				{
					apply();
					compute();
					if(equilibrium() || !corrector()) break;
				}
				if(m_status) break;
				restore();
				compute();
			}
			if(m_energy) m_energy(m_U, m_V, m_x_new, m_p_new);
			update();
			record();
		}

		//data
		uint32_t NewtonRaphson::state_set(void) const
		{
			return 1 << uint32_t(State::x) | 1 << uint32_t(State::p) | 1 << uint32_t(State::t);
		}
		uint32_t NewtonRaphson::force_set(void) const
		{
			return 1 << uint32_t(Force::r) | 1 << uint32_t(Force::fe);
		}
		uint32_t NewtonRaphson::energy_set(void) const
		{
			return 1 << uint32_t(Energy::U) | 1 << uint32_t(Energy::V);
		}
		uint32_t NewtonRaphson::tangent_set(void) const
		{
			return 1 << uint32_t(Tangent::K);
		}

		//data
		NewtonRaphson::SystemCallback NewtonRaphson::system(void) const
		{
			return m_system;
		}
		NewtonRaphson::SystemCallback NewtonRaphson::system(NewtonRaphson::SystemCallback system)
		{
			return m_system = system;
		}

		NewtonRaphson::EnergyCallback NewtonRaphson::energy(void) const
		{
			return m_energy;
		}
		NewtonRaphson::EnergyCallback NewtonRaphson::energy(NewtonRaphson::EnergyCallback energy)
		{
			return m_energy = energy;
		}

		//solve
		void NewtonRaphson::check(void)
		{
			if(!m_system)
			{
				throw std::runtime_error("Newton-Raphson solver called with system unset!");
			}
		}
		void NewtonRaphson::print(void)
		{
			Incremental::print();
			Implicit::print();
			Solver::print();
		}
		void NewtonRaphson::setup(void)
		{
			Solver::setup();
			Implicit::setup();
			Incremental::setup();
			if(m_energy) m_energy(m_U, m_V, m_x_new, m_p_new);
		}
		void NewtonRaphson::compute(void)
		{
			m_system(m_r, m_fe, m_K, m_x_new, m_p_new);
		}
		bool NewtonRaphson::predictor(void)
		{
			//predictor
			m_dp = m_dps;
			memcpy(m_dx, m_dxs, m_size * sizeof(double));
			if(!Solver::solve(m_K, m_r, m_dxr) || !Solver::solve(m_K, m_fe, m_dxt))
			{
				if(!m_silent) printf("Unable to decompose stiffness matrix in predictor!\n");
				return false;
			}
			//continuation
			m_dp = (m_step == 1? m_dps : m_continuation.predictor()) / (1 << m_attempt);
			for(uint32_t i = 0; i < m_size; i++) m_dx[i] = m_ddx[i] = m_dxr[i] + m_dp * m_dxt[i];
			//return
			return true;
		}
		bool NewtonRaphson::corrector(void)
		{
			//corrector
			if(!Solver::solve(m_K, m_r, m_ddxr) || !Solver::solve(m_K, m_fe, m_ddxt))
			{
				if(!m_silent) printf("Unable to decompose stiffness matrix in corrector!\n");
				return false;
			}
			m_ddp = m_continuation.corrector();
			//update
			m_dp += m_ddp;
			for(uint32_t i = 0; i < m_size; i++) m_dx[i] += m_ddx[i] = m_ddxr[i] + m_ddp * m_ddxt[i];
			//return
			return true;
		}
	}
}