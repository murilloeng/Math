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
			//update
			update();
			//record
			record();
		}

		//data
		uint32_t RungeKutta::state_set(void) const
		{
			return 1 << (uint32_t) State::x | 1 << (uint32_t) State::v | 1 << (uint32_t) State::a | 1 << (uint32_t) State::t;
		}
		uint32_t RungeKutta::force_set(void) const
		{
			return 1 << (uint32_t) Force::fi | 1 << (uint32_t) Force::fe;
		}
		uint32_t RungeKutta::tangent_set(void) const
		{
			return 1 << (uint32_t) Tangent::M;
		}

		//data
		RungeKutta::Inertia RungeKutta::inertia(void) const
		{
			return m_inertia;
		}
		RungeKutta::Inertia RungeKutta::inertia(Inertia inertia)
		{
			return m_inertia = inertia;
		}

		RungeKutta::InternalForce RungeKutta::internal_force(void) const
		{
			return m_internal_force;
		}
		RungeKutta::InternalForce RungeKutta::internal_force(InternalForce internal_force)
		{
			return m_internal_force = internal_force;
		}

		RungeKutta::ExternalForce RungeKutta::external_force(void) const
		{
			return m_external_force;
		}
		RungeKutta::ExternalForce RungeKutta::external_force(ExternalForce external_force)
		{
			return m_external_force = external_force;
		}

		//solve
		void RungeKutta::check(void)
		{
			if(!m_inertia || !m_internal_force || !m_external_force)
			{
				throw std::runtime_error("Runge Kutta solver called with at least one method not set!");
			}
		}
		void RungeKutta::print(void)
		{
			Incremental::print();
			Solver::print();
		}
		void RungeKutta::compute(void)
		{
			//setup
			m_inertia(m_M, m_x_new);
			m_internal_force(m_fi, m_x_new, m_v_new);
			m_external_force(m_fe, m_x_new, m_v_new, m_t_new);
			for(uint32_t i = 0; i < m_size; i++) m_fe[i] -= m_fi[i];
			//compute
			Solver::solve(m_M, m_fe, m_a_new);
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