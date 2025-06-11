//std
#include <cstdio>
#include <cstring>

//math
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/newmark.hpp"

namespace math
{
	namespace solvers
	{
		//constructors
		newmark::newmark(void) : 
			m_args(nullptr), 
			m_watch_dof(0), m_size(0), 
			m_step_max(0), m_attempt_max(5), m_iteration_max(10), 
			m_T(1), m_g(0.50), m_b(0.25),
			m_x_old(nullptr), m_x_new(nullptr), m_x_data(nullptr), m_dx(nullptr),
			m_v_old(nullptr), m_v_new(nullptr), m_v_data(nullptr), m_dv(nullptr),
			m_a_old(nullptr), m_a_new(nullptr), m_a_data(nullptr), m_da(nullptr),
			m_r(nullptr), m_fe(nullptr), m_fi(nullptr), m_M(nullptr), m_C(nullptr), m_K(nullptr),
			m_internal_force(nullptr), m_external_force(nullptr),
			m_inertia(nullptr), m_damping(nullptr), m_stiffness(nullptr)
		{
			return;
		}

		//destructor
		newmark::~newmark(void)
		{
			cleanup();
		}

		//data
		// void newmark::update(void)
		// {
		// 	// inertia();
		// 	// internal();
		// 	// external();
		// 	// for(uint32_t i = 0; i < m_nd; i++)
		// 	// {
		// 	// 	m_r[i] = m_fe[i] - m_fi[i];
		// 	// }
		// 	// math::vector am(m_a, m_nd);
		// 	// math::vector rm(m_r, m_nd);
		// 	// math::matrix(m_M, m_nd, m_nd).solve(am, rm);
		// }
		void newmark::residue(void)
		{
			// inertia();
			// internal();
			// external();
			// for(uint32_t i = 0; i < m_nd; i++)
			// {
			// 	m_r[i] = m_fe[i] - m_fi[i];
			// 	for(uint32_t j = 0; j < m_nd; j++)
			// 	{
			// 		m_r[i] -= m_M[i + m_nd * j] * m_a[j];
			// 	}
			// }
		}
		void newmark::inertia(void)
		{
			m_inertia(m_M, m_x_new, m_args);
		}
		void newmark::damping(void)
		{
			m_damping(m_C, m_x_new, m_v_new, m_t, m_args);
		}
		void newmark::stifness(void)
		{
			m_stiffness(m_K, m_x_new, m_v_new, m_a_new, m_t, m_args);
		}
		void newmark::internal(void)
		{
			m_internal_force(m_fi, m_x_new, m_v_new, m_args);
		}
		void newmark::external(void)
		{
			m_external_force(m_fe, m_x_new, m_v_new, m_t, m_args);
		}

		//solve
		bool newmark::stop(void)
		{
			return false;
		}
		void newmark::check(void)
		{

		}
		void newmark::apply(void)
		{

		}
		void newmark::print(void)
		{

		}
		void newmark::setup(void)
		{
			m_t = 0;
			m_step = 0;
			m_dt = m_T / m_step_max;
		}
		void newmark::predictor(void)
		{
			m_t += m_dt;
			for(uint32_t i = 0; i < m_size; i++)
			{
				m_v_new[i] = m_v_old[i] + m_dt * m_a_old[i];
				m_x_new[i] = m_x_old[i] + m_dt * m_v_old[i] - m_dt * m_dt / 2 * m_a_old[i];
			}
		}
		void newmark::corrector(void)
		{
			// while(true)
			// {
			// 	residue();
			// 	math::vector rm(m_r, m_size), dxm(m_dx, m_size);
			// 	const double f = math::vector(m_fi, m_size).norm();
			// 	if(math::vector(m_r, m_size).norm() < 1e-5 * (f == 0 ? 1 : f))
			// 	{
			// 		break;
			// 	}
			// 	damping();
			// 	stifness();
			// 	for(uint32_t i = 0; i < m_size * m_size; i++)
			// 	{
			// 		m_K[i] += (m_g * m_dt * m_C[i] + m_M[i]) / (m_b * m_dt * m_dt);
			// 	}
			// 	math::matrix(m_K, m_size, m_size).solve(dxm, rm);
			// 	for(uint32_t i = 0; i < m_size; i++)
			// 	{
			// 		m_x[i] += m_dx[i];
			// 		m_v[i] += m_dx[i] * m_g / (m_b * m_dt);
			// 		m_a[i] += m_dx[i] / (m_b * m_dt * m_dt);
			// 	}
			// }
		}
		bool newmark::equilibrium(void)
		{
			return false;
		}

		//solve
		void newmark::step(void)
		{
			// predictor();
			// corrector();
		}
		void newmark::solve(void)
		{
			// setup();
			// update();
			// serialize();
			// while(m_step < m_step_max)
			// {
			// 	step();
			// 	serialize();
			// }
		}
		void newmark::cleanup(void)
		{
			//data
			double** data[] = {
				&m_x_old, &m_x_new, &m_x_data, &m_dx,
				&m_v_old, &m_v_new, &m_v_data, &m_dv,
				&m_a_old, &m_a_new, &m_a_data, &m_da,
				&m_r, &m_fe, &m_fi, &m_M, &m_C, &m_K
			};
			//delete
			for(double** ptr : data)
			{
				delete[] *ptr;
				*ptr = nullptr;
			}
		}
		void newmark::allocate(void)
		{
			m_r = new double[m_size];
			m_fe = new double[m_size];
			m_fi = new double[m_size];
			m_dx = new double[m_size];
			m_dv = new double[m_size];
			m_da = new double[m_size];
			m_x_old = new double[m_size];
			m_x_new = new double[m_size];
			m_v_old = new double[m_size];
			m_v_new = new double[m_size];
			m_a_old = new double[m_size];
			m_a_new = new double[m_size];
			m_M = new double[m_size * m_size];
			m_C = new double[m_size * m_size];
			m_K = new double[m_size * m_size];
			m_x_data = new double[m_size * (m_step_max + 1)];
			m_v_data = new double[m_size * (m_step_max + 1)];
			m_a_data = new double[m_size * (m_step_max + 1)];
		}
	}
}