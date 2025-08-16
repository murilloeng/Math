//std
#include <cstdio>
#include <cstring>
#include <stdexcept>

//math
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/newmark.hpp"

namespace math
{
	namespace solvers
	{
		//constructors
		newmark::newmark(void) : m_g(0.50), m_b(0.25)
		{
			return;
		}

		//destructor
		newmark::~newmark(void)
		{
			cleanup();
		}

		//data
		uint32_t newmark::state_set(void) const
		{
			return uint32_t(state::x) | uint32_t(state::v) | uint32_t(state::a) | uint32_t(state::t);
		}
		uint32_t newmark::force_set(void) const
		{
			return uint32_t(force::r) | uint32_t(force::fi) | uint32_t(force::fe);
		}
		uint32_t newmark::tangent_set(void) const
		{
			return uint32_t(tangent::K) | uint32_t(tangent::C) | uint32_t(tangent::M);
		}

		//solve
		void newmark::check(void)
		{
			if(!m_internal_force || !m_external_force || !m_inertia || !m_damping || !m_stiffness)
			{
				throw std::runtime_error("Newmark solver called with at least one method not set!");
			}
		}
		void newmark::setup(void)
		{
			//data
			vector a(m_a_new, m_size);
			const vector r(m_r, m_size);
			const matrix M(m_M, m_size, m_size);
			//setup
			solver::setup();
			for(uint32_t i = 0; i < m_size; i++)
			{
				m_r[i] = m_fe[i] - m_fi[i];
			}
			//check
			if(!M.solve(a, r))
			{
				if(!m_silent) printf("Unable to compute acceleration in setup!\n");
			}
			memcpy(m_a_old, m_a_new, m_size * sizeof(double));
		}
		void newmark::compute(void)
		{
			//forces
			m_internal_force(m_fi, m_x_new, m_v_new);
			m_external_force(m_fe, m_x_new, m_v_new, m_t_new);
			//tangents
			m_inertia(m_M, m_x_new);
			m_damping(m_C, m_x_new, m_v_new, m_t_new);
			m_stiffness(m_K, m_x_new, m_v_new, m_a_new, m_t_new);
			//residue
			for(uint32_t i = 0; i < m_size; i++)
			{
				m_r[i] = m_fe[i] - m_fi[i];
				for(uint32_t j = 0; j < m_size; j++)
				{
					m_r[i] -= m_M[i + m_size * j] * m_a_new[j];
				}
			}
		}
		void newmark::predictor(void)
		{
			for(uint32_t i = 0; i < m_size; i++)
			{
				m_da[i] = 0;
				m_dv[i] = m_dt * m_a_old[i];
				m_dx[i] = m_dt * m_v_old[i] + m_dt * m_dt / 2 * m_a_old[i];
			}
		}
		void newmark::corrector(void)
		{
			//data
			vector ddxr(m_ddxr, m_size);
			const vector r(m_r, m_size);
			const matrix K(m_K, m_size, m_size);
			//tangent
			for(uint32_t i = 0; i < m_size * m_size; i++)
			{
				m_K[i] += m_g * m_C[i] / m_b / m_dt + m_M[i] / m_b / m_dt / m_dt;
			}
			//solve
			if(!K.solve(ddxr, r))
			{
				if(!m_silent) printf("Unable to decompose system matrix in corrector!\n");
			}
			//update
			for(uint32_t i = 0; i < m_size; i++)
			{
				m_dx[i] += m_ddxr[i];
				m_dv[i] += m_g * m_ddxr[i] / m_b / m_dt;
				m_da[i] += m_ddxr[i] / m_b / m_dt / m_dt;
			}
		}
	}
}