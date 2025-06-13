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
				printf("Error: Newmark solver called with methods not set!\n");
				exit(EXIT_FAILURE);
			}
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
			apply();
		}
		void newmark::corrector(void)
		{
			//v_new = v_old + dt * a_old + g * dt * (a_new - a_old)
			//x_new = x_old + dt * v_old + dt * dt / 2 * a_old + b * dt * dt * (a_new - a_old)

			//r(x(t), v(t), a(t)) = fe(x(t), v(t), t) - fi(x(t), v(t)) - M(x(t)) * a(t) = 0

			//dv_new = g * dt * da_new
			//dx_new = b * dt * dt * da_new

			//r - K * dx - C * dv - M * da = 0
			//r - (K + g / b / dt * C + 1 / b / dt / dt * M) * dx = 0

			//S = K + g / b / dt * C + 1 / b / dt / dt * M

			//dx = S^{-1} * r
			//dv = g / b / dt * dx
			//da = 1 / b / dt / dt * dx
		}
	}
}