//std
#include <cstdio>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Solvers/Newmark.hpp"

//x: state Vector
//r: residue Vector
//v: velocity Vector
//a: acceleration Vector
//p: continuation parameter

//target system: r(x(t), v(t), a(t)) = fe(x(t), v(t), t) - fi(x(t), v(t)) - M(x(t)) * a(t) = 0

//tangent on a: M(x) = -dr/da(x)
//tangent on v: C(x, v, t) = -dr/dv(x, v, t) = dfi/dv(x, v) - dfe/dv(x, v, t)
//tangent on x: K(x, v, a) = -dr/dx(x, v, a, t) = dfi/dx(x, v) - dfe/dx(x, v, t) + dM/dx(x) : a

//time interpolation:
//v_new = v_old + dt * a_old + g * dt * (a_new - a_old)
//x_new = x_old + dt * v_old + dt * dt / 2 * a_old + b * dt * dt * (a_new - a_old)

namespace math
{
	namespace solvers
	{
		//constructors
		Newmark::Newmark(void) : m_g(0.50), m_b(0.25)
		{
			return;
		}

		//destructor
		Newmark::~Newmark(void)
		{
			return;
		}

		//data
		uint32_t Newmark::state_set(void) const
		{
			return uint32_t(State::x) | uint32_t(State::v) | uint32_t(State::a) | uint32_t(State::t);
		}
		uint32_t Newmark::force_set(void) const
		{
			return uint32_t(Force::r) | uint32_t(Force::fi) | uint32_t(Force::fe);
		}
		uint32_t Newmark::tangent_set(void) const
		{
			return uint32_t(Tangent::K) | uint32_t(Tangent::C) | uint32_t(Tangent::M);
		}

		//solve
		void Newmark::check(void)
		{
			if(!m_internal_force || !m_external_force || !m_inertia || !m_damping || !m_stiffness)
			{
				throw std::runtime_error("Newmark solver called with at least one method not set!");
			}
		}
		void Newmark::setup(void)
		{
			//setup
			Solver::setup();
			for(uint32_t i = 0; i < m_size; i++)
			{
				m_r[i] = m_fe[i] - m_fi[i];
			}
			//check
			if(!solve(m_M, m_r, m_a_new))
			{
				if(!m_silent) printf("Unable to compute acceleration in setup!\n");
			}
			memcpy(m_a_old, m_a_new, m_size * sizeof(double));
		}
		void Newmark::compute(void)
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
		void Newmark::predictor(void)
		{
			for(uint32_t i = 0; i < m_size; i++)
			{
				m_da[i] = 0;
				m_dv[i] = m_dt * m_a_old[i];
				m_dx[i] = m_dt * m_v_old[i] + m_dt * m_dt / 2 * m_a_old[i];
			}
		}
		void Newmark::corrector(void)
		{
			//tangent
			for(uint32_t i = 0; i < m_size * m_size; i++)
			{
				m_K[i] += m_g * m_C[i] / m_b / m_dt + m_M[i] / m_b / m_dt / m_dt;
			}
			//solve
			if(!solve(m_K, m_r, m_ddxr))
			{
				if(!m_silent) printf("Unable to decompose system Matrix in corrector!\n");
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