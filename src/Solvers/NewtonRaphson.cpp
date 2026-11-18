//std
#include <cmath>
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
			if(m_step_adapt && m_step > 1) compute_step_scale();
			compute_step_clip();
			compute_step_solve();
			if(!m_status) return;
			if(m_step_adapt) m_step_scale /= 1 << m_attempt;
			compute_step_limit();
			if(m_step_stop) return;
			if(m_energy) m_energy(m_U, m_V, m_x_new, m_p_new);
			update();
			record();
		}
		void NewtonRaphson::compute_step_clip(void)
		{
			//data
			const uint32_t types = m_stop_criteria.types();
			const double p_min = m_stop_criteria.load_min();
			const double p_max = m_stop_criteria.load_max();
			const bool test_min = types & 1 << uint32_t(StopCriteria::Type::LoadLimitMinimum);
			const bool test_max = types & 1 << uint32_t(StopCriteria::Type::LoadLimitMaximum);
			//clip
			if(m_step == 1 || m_continuation.type() == Continuation::Type::LoadControl)
			{
				if(test_min && m_p_old + m_dps < p_min) m_dps = p_min - m_p_old;
				if(test_max && m_p_old + m_dps > p_max) m_dps = p_max - m_p_old;
			}
		}
		void NewtonRaphson::compute_step_solve(void)
		{
			for(m_attempt = 0; m_attempt < m_attempt_max; m_attempt++)
			{
				// if(m_step_adapt && m_step_scale / (1 << m_attempt) < m_step_scale_min) break;
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
		}
		void NewtonRaphson::compute_step_limit(void)
		{
			//data
			const uint32_t types = m_stop_criteria.types();
			const double p_min = m_stop_criteria.load_min();
			const double p_max = m_stop_criteria.load_max();
			const bool test_min = types & 1 << uint32_t(StopCriteria::Type::LoadLimitMinimum);
			const bool test_max = types & 1 << uint32_t(StopCriteria::Type::LoadLimitMaximum);
			const double t_min = 1e-12 * fmax(1, fabs(p_min));
			const double t_max = 1e-12 * fmax(1, fabs(p_max));
			//landed
			if(test_min && fabs(m_p_new - p_min) <= t_min) { m_p_new = p_min; return; }
			if(test_max && fabs(m_p_new - p_max) <= t_max) { m_p_new = p_max; return; }
			//check
			const bool over_min = test_min && m_p_new < p_min;
			const bool over_max = test_max && m_p_new > p_max;
			if(!over_min && !over_max) return;
			//stop (step not recorded)
			restore();
			compute();
			m_step_stop = true;
		}
		void NewtonRaphson::compute_step_scale(void)
		{
			//factor
			const double r = double(m_iteration_desired) / (m_iteration > 1 ? m_iteration : 1);
			double f = fmax(m_step_factor_min, fmin(m_step_factor_max, sqrt(r)));
			if(m_attempt != 0) f = fmin(f, 1);
			//bounds
			f = fmax(m_step_scale_min, fmin(m_step_scale_max, m_step_scale * f)) / m_step_scale;
			//scale
			m_dps *= f;
			m_step_scale *= f;
			for(uint32_t i = 0; i < m_size; i++) m_dxs[i] *= f;
		}

		//data
		uint32_t NewtonRaphson::states_set(void) const
		{
			return 1 << uint32_t(State::x) | 1 << uint32_t(State::p) | 1 << uint32_t(State::t);
		}
		uint32_t NewtonRaphson::forces_set(void) const
		{
			return 1 << uint32_t(Force::r) | 1 << uint32_t(Force::fe);
		}
		uint32_t NewtonRaphson::energies_set(void) const
		{
			return 1 << uint32_t(Energy::U) | 1 << uint32_t(Energy::V);
		}
		uint32_t NewtonRaphson::tangents_set(void) const
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