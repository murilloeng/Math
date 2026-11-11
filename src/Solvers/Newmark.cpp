//std
#include <cstdio>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Solvers/Newmark.hpp"

//x: state
//r: residue
//v: velocity
//a: acceleration
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
		Newmark::Newmark(void) : m_g{0.50}, m_b{0.25}
		{
			return;
		}

		//solve
		void Newmark::compute_step(void)
		{
			m_status = false;
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
			}
			if(m_energy) m_energy(m_T, m_U, m_V, m_x_new, m_v_new, m_t_new);
			update();
			record();
		}

		//data
		uint32_t Newmark::states_set(void) const
		{
			return 1 << uint32_t(State::x) | 1 << uint32_t(State::v) | 1 << uint32_t(State::a) | 1 << uint32_t(State::t);
		}
		uint32_t Newmark::energies_set(void) const
		{
			return 1 << uint32_t(Energy::T) | 1 << uint32_t(Energy::U) | 1 << uint32_t(Energy::V);
		}
		uint32_t Newmark::tangents_set(void) const
		{
			return 1 << uint32_t(Tangent::K) | 1 << uint32_t(Tangent::C) | 1 << uint32_t(Tangent::M);
		}

		//data
		double Newmark::beta(void) const
		{
			return m_b;
		}
		double Newmark::beta(double beta)
		{
			return m_b = beta;
		}

		double Newmark::gamma(void) const
		{
			return m_g;
		}
		double Newmark::gamma(double gamma)
		{
			return m_g = gamma;
		}

		Newmark::EnergyCallback Newmark::energy(void) const
		{
			return m_energy;
		}
		Newmark::EnergyCallback Newmark::energy(EnergyCallback energy)
		{
			return m_energy = energy;
		}

		Newmark::SystemCallback Newmark::system(void) const
		{
			return m_system;
		}
		Newmark::SystemCallback Newmark::system(SystemCallback system)
		{
			return m_system = system;
		}

		//solve
		void Newmark::check(void)
		{
			if(!m_system)
			{
				throw std::runtime_error("Newmark solver called with system unset!");
			}
		}
		void Newmark::print(void)
		{
			Incremental::print();
			Implicit::print();
			Solver::print();
		}
		void Newmark::setup(void)
		{
			//setup
			Solver::setup();
			Implicit::setup();
			Incremental::setup();
			m_system(m_fi, m_fe, m_K, m_C, m_M, m_x_new, m_v_new, m_a_new, m_t_new);
			//residue
			for(uint32_t i = 0; i < m_size; i++)
			{
				m_r[i] = m_fe[i] - m_fi[i];
			}
			//acceleration
			if(!Solver::solve(m_M, m_r, m_a_new))
			{
				if(!m_silent) printf("Unable to compute acceleration in setup!\n");
			}
			memcpy(m_a_old, m_a_new, m_size * sizeof(double));
			//energy
			if(m_energy) m_energy(m_T, m_U, m_V, m_x_new, m_v_new, m_t_new);
		}
		void Newmark::compute(void)
		{
			//compute
			m_system(m_fi, m_fe, m_K, m_C, m_M, m_x_new, m_v_new, m_a_new, m_t_new);
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
		bool Newmark::predictor(void)
		{
			//predictor
			for(uint32_t i = 0; i < m_size; i++)
			{
				m_da[i] = 0;
				m_dv[i] = m_dt * m_a_old[i];
				m_dx[i] = m_ddx[i] = m_dt * m_v_old[i] + m_dt * m_dt / 2 * m_a_old[i];
			}
			//return
			return true;
		}
		bool Newmark::corrector(void)
		{
			//tangent
			for(uint32_t i = 0; i < m_size * m_size; i++)
			{
				m_K[i] += m_g * m_C[i] / m_b / m_dt + m_M[i] / m_b / m_dt / m_dt;
			}
			//solve
			if(!Solver::solve(m_K, m_r, m_ddx))
			{
				if(!m_silent) printf("Unable to decompose system Matrix in corrector!\n");
				return false;
			}
			//update
			for(uint32_t i = 0; i < m_size; i++)
			{
				m_dx[i] += m_ddx[i];
				m_dv[i] += m_g * m_ddx[i] / m_b / m_dt;
				m_da[i] += m_ddx[i] / m_b / m_dt / m_dt;
			}
			//return
			return true;
		}
	}
}