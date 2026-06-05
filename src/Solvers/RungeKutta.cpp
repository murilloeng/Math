//std
#include <cstdio>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Linear/Vector.hpp"
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
			//data
			math::Vector a(m_a_new, m_size);
			math::Matrix M(m_M, m_size, m_size);
			math::Vector fi(m_fi, m_size), fe(m_fe, m_size);
			//setup
			m_inertia(m_M, m_x_new);
			m_internal_force(m_fi, m_x_new, m_v_new);
			m_external_force(m_fe, m_x_new, m_v_new, m_t_new);
			//compute
			M.solve(a, fe - fi);
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

		// //data
		// void RungeKutta::update(void)
		// {
		// 	!m_type ? m_system_1(m_v, m_x, m_t) : m_system_2(m_a, m_x, m_v, m_t);
		// }
		// void RungeKutta::state(double f)
		// {
		// 	if(!m_type)
		// 	{
		// 		for(uint32_t i = 0; i < m_nd; i++)
		// 		{
		// 			m_x[i] = m_xn[i] + f * m_v[i];
		// 		}
		// 	}
		// 	else
		// 	{
		// 		for(uint32_t i = 0; i < m_nd; i++)
		// 		{
		// 			m_x[i] = m_xn[i] + f * m_v[i];
		// 			m_v[i] = m_vn[i] + f * m_a[i];
		// 		}
		// 	}
		// }
		// void RungeKutta::increment(double f)
		// {
		// 	if(!m_type)
		// 	{
		// 		for(uint32_t i = 0; i < m_nd; i++)
		// 		{
		// 			m_dx[i] += f * m_v[i];
		// 		}
		// 	}
		// 	else
		// 	{
		// 		for(uint32_t i = 0; i < m_nd; i++)
		// 		{
		// 			m_dx[i] += f * m_v[i];
		// 			m_dv[i] += f * m_a[i];
		// 		}
		// 	}
		// }

		// //solve
		// void RungeKutta::setup(void)
		// {
		// 	m_t = 0;
		// 	m_s = 0;
		// 	m_dt = m_T / m_ns;
		// 	memcpy(m_xn, m_x, m_nd * sizeof(double));
		// 	memcpy(m_vn, m_v, m_nd * sizeof(double));
		// }
		// void RungeKutta::tangent_1(void)
		// {
		// 	memset(m_dx, 0, m_nd * sizeof(double));
		// 	memset(m_dv, 0, m_nd * sizeof(double));
		// 	increment(m_dt / 6);
		// }
		// void RungeKutta::tangent_2(void)
		// {
		// 	//time
		// 	m_t += m_dt / 2;
		// 	//update state
		// 	state(m_dt / 2);
		// 	//update tangent
		// 	update();
		// 	increment(m_dt / 3);
		// }
		// void RungeKutta::tangent_3(void)
		// {
		// 	//update state
		// 	state(m_dt / 2);
		// 	//update tangent
		// 	update();
		// 	increment(m_dt / 3);
		// }
		// void RungeKutta::tangent_4(void)
		// {
		// 	//time
		// 	m_t += m_dt / 2;
		// 	//update state
		// 	state(m_dt);
		// 	//update tangent
		// 	update();
		// 	increment(m_dt / 6);
		// }
		// void RungeKutta::corrector(void)
		// {
		// 	//update state
		// 	for(uint32_t i = 0; i < m_nd; i++)
		// 	{
		// 		m_x[i] = m_xn[i] += m_dx[i];
		// 		m_v[i] = m_vn[i] += m_dv[i];
		// 	}
		// 	//update system
		// 	update();
		// }

		// //solve
		// void RungeKutta::init(void)
		// {
		// 	setup();
		// 	update();
		// 	serialize();
		// }
		// void RungeKutta::step(void)
		// {
		// 	tangent_1();
		// 	tangent_2();
		// 	tangent_3();
		// 	tangent_4();
		// 	corrector();
		// }
		// void RungeKutta::solve(void)
		// {
		// 	init();
		// 	while(m_s < m_ns)
		// 	{
		// 		step();
		// 		serialize();
		// 	}
		// }
		// void RungeKutta::serialize(void)
		// {
		// 	m_s++;
		// 	printf("%04d %+.6e ", m_s, m_t);
		// 	for(uint32_t i = 0; i < m_nd; i++)
		// 	{
		// 		printf("%+.6e %+.6e %+.6e ", m_x[i], m_v[i], m_a[i]);
		// 	}
		// 	printf("\n");
		// }
	}
}