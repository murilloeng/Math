//std
#include <cstdio>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Solvers/RungeKutta.hpp"

namespace math
{
	namespace solvers
	{
		//solve
		void RungeKutta::compute_step(void)
		{
			//predictor
			predictor();
			//corrector
			apply();
			compute();
			corrector();
			m_status = true;
			if(m_energy) m_energy(m_T, m_U, m_V, m_x_new, m_v_new, m_t_new);
			//update
			update();
			//record
			record();
		}

		//data
		uint32_t RungeKutta::state_set(void) const
		{
			return 1 << uint32_t(State::x) | 1 << uint32_t(State::v) | 1 << uint32_t(State::a) | 1 << uint32_t(State::t);
		}
		uint32_t RungeKutta::force_set(void) const
		{
			return 1 << uint32_t(Force::fi) | 1 << uint32_t(Force::fe);
		}
		uint32_t RungeKutta::energy_set(void) const
		{
			return 1 << uint32_t(Energy::T) | 1 << uint32_t(Energy::U) | 1 << uint32_t(Energy::V);
		}
		uint32_t RungeKutta::tangent_set(void) const
		{
			return 1 << uint32_t(Tangent::M);
		}

		//data
		RungeKutta::EnergyCallback RungeKutta::energy(void) const
		{
			return m_energy;
		}
		RungeKutta::EnergyCallback RungeKutta::energy(EnergyCallback energy)
		{
			return m_energy = energy;
		}

		RungeKutta::SystemCallback RungeKutta::system(void) const
		{
			return m_system;
		}
		RungeKutta::SystemCallback RungeKutta::system(SystemCallback system)
		{
			return m_system = system;
		}

		//solve
		void RungeKutta::check(void)
		{
			if(!m_system)
			{
				throw std::runtime_error("Runge Kutta solver called with system unset!");
			}
		}
		void RungeKutta::print(void)
		{
			Incremental::print();
			Solver::print();
		}
		void RungeKutta::setup(void)
		{
			//setup
			Solver::setup();
			Incremental::setup();
			m_system(m_M, m_fi, m_fe, m_x_new, m_v_new, m_t_new);
			//residue
			for(uint32_t i = 0; i < m_size; i++)
			{
				m_fe[i] -= m_fi[i];
			}
			//acceleration
			if(!Solver::solve(m_M, m_fe, m_a_new))
			{
				if(!m_silent) printf("Unable to compute acceleration in setup!\n");
			}
			memcpy(m_a_old, m_a_new, m_size * sizeof(double));
			//energy
			if(m_energy) m_energy(m_T, m_U, m_V, m_x_new, m_v_new, m_t_new);
		}
		void RungeKutta::compute(void)
		{
			//setup
			m_system(m_M, m_fi, m_fe, m_x_new, m_v_new, m_t_new);
			for(uint32_t i = 0; i < m_size; i++) m_fe[i] -= m_fi[i];
			//compute
			Solver::solve(m_M, m_fe, m_a_new);
		}
		bool RungeKutta::predictor(void)
		{
			memset(m_dx, 0, m_size * sizeof(double));
			memset(m_dv, 0, m_size * sizeof(double));
			return true;
		}
		bool RungeKutta::corrector(void)
		{
			compute_tangent_1();
			compute_tangent_2();
			compute_tangent_3();
			compute_tangent_4();
			return true;
		}

		//compute
		void RungeKutta::compute_tangent_1(void)
		{
			//setup
			m_t_new = m_t_old;
			for(uint32_t i = 0; i < m_size; i++) m_x_new[i] = m_x_old[i];
			for(uint32_t i = 0; i < m_size; i++) m_v_new[i] = m_v_old[i];
			//update
			compute();
			for(uint32_t i = 0; i < m_size; i++) m_dx[i] += m_dt / 6 * m_v_new[i];
			for(uint32_t i = 0; i < m_size; i++) m_dv[i] += m_dt / 6 * m_a_new[i];
		}
		void RungeKutta::compute_tangent_2(void)
		{
			//setup
			m_t_new = m_t_old + m_dt / 2;
			for(uint32_t i = 0; i < m_size; i++) m_x_new[i] = m_x_old[i] + m_dt / 2 * m_v_new[i];
			for(uint32_t i = 0; i < m_size; i++) m_v_new[i] = m_v_old[i] + m_dt / 2 * m_a_new[i];
			//update
			compute();
			for(uint32_t i = 0; i < m_size; i++) m_dx[i] += m_dt / 3 * m_v_new[i];
			for(uint32_t i = 0; i < m_size; i++) m_dv[i] += m_dt / 3 * m_a_new[i];
		}
		void RungeKutta::compute_tangent_3(void)
		{
			//setup
			m_t_new = m_t_old + m_dt / 2;
			for(uint32_t i = 0; i < m_size; i++) m_x_new[i] = m_x_old[i] + m_dt / 2 * m_v_new[i];
			for(uint32_t i = 0; i < m_size; i++) m_v_new[i] = m_v_old[i] + m_dt / 2 * m_a_new[i];
			//update
			compute();
			for(uint32_t i = 0; i < m_size; i++) m_dx[i] += m_dt / 3 * m_v_new[i];
			for(uint32_t i = 0; i < m_size; i++) m_dv[i] += m_dt / 3 * m_a_new[i];
		}
		void RungeKutta::compute_tangent_4(void)
		{
			//setup
			m_t_new = m_t_old + m_dt;
			for(uint32_t i = 0; i < m_size; i++) m_x_new[i] = m_x_old[i] + m_dt * m_v_new[i];
			for(uint32_t i = 0; i < m_size; i++) m_v_new[i] = m_v_old[i] + m_dt * m_a_new[i];
			//update
			compute();
			for(uint32_t i = 0; i < m_size; i++) m_dx[i] += m_dt / 6 * m_v_new[i];
			for(uint32_t i = 0; i < m_size; i++) m_dv[i] += m_dt / 6 * m_a_new[i];
		}
	}
}