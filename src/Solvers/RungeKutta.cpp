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
		//constructors
		RungeKutta::RungeKutta(void)
		{
			return;
		}

		//destructor
		RungeKutta::~RungeKutta(void)
		{
			return;
		}

		//data
		uint32_t RungeKutta::state_set(void) const
		{
			return 
				(uint32_t) State::x |
				(uint32_t) State::v |
				(uint32_t) State::a |
				(uint32_t) State::t;
		}
		uint32_t RungeKutta::force_set(void) const
		{
			return
				(uint32_t) Force::fi |
				(uint32_t) Force::fe;
		}
		uint32_t RungeKutta::tangent_set(void) const
		{
			return (uint32_t) Tangent::M;
		}

		//solve
		void RungeKutta::check(void)
		{
			if(!m_inertia || !m_internal_force || !m_external_force)
			{
				throw std::runtime_error("Runge Kutta solver called with at least one method not set!");
			}
		}
		void RungeKutta::compute(void)
		{
			//setup
			m_inertia(m_M, m_x_new);
			m_internal_force(m_fi, m_x_new, m_v_new);
			m_external_force(m_fe, m_x_new, m_v_new, m_t_new);
			for(uint32_t i = 0; i < m_size; i++) m_fe[i] -= m_fi[i];
			//compute
			solve(m_M, m_fe, m_a_new);
		}
		void RungeKutta::predictor(void)
		{
			memset(m_dx, 0, m_size * sizeof(double));
			memset(m_dv, 0, m_size * sizeof(double));
		}
		void RungeKutta::corrector(void)
		{
			compute_tangent_1();
			compute_tangent_2();
			compute_tangent_3();
			compute_tangent_4();
		}
		bool RungeKutta::equilibrium(void)
		{
			return m_equilibrium = m_iteration == 1;
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