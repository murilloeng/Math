//std
#include <cstdio>
#include <stdexcept>

//Math
#include "Math/inc/Solvers/NewtonRaphson.hpp"

//x: state Vector
//r: residual force Vector
//p: continuation parameter

//target system: r(x, p) = 0

//tangent on x: K(x, p) = -dr/dx(x, p)
//tangent on p: g(x, p) = +dr/dp(x, p)

namespace math
{
	namespace solvers
	{
		//constructors
		NewtonRaphson::NewtonRaphson(void)
		{
			return;
		}

		//destructor
		NewtonRaphson::~NewtonRaphson(void)
		{
			return;
		}

		//data
		uint32_t NewtonRaphson::state_set(void) const
		{
			return uint32_t(State::x) | uint32_t(State::p);
		}
		uint32_t NewtonRaphson::force_set(void) const
		{
			return uint32_t(Force::r) | uint32_t(Force::fe);
		}
		uint32_t NewtonRaphson::tangent_set(void) const
		{
			return uint32_t(Tangent::K);
		}

		//solve
		void NewtonRaphson::check(void)
		{
			if(!m_system_1 && !m_system_2 && !(m_residue && m_tangent_1 && m_tangent_2))
			{
				throw std::runtime_error("Newton-Raphson solver called with at least one method not set!");
			}
		}
		void NewtonRaphson::compute(void)
		{
			if(m_system_2)
			{
				m_system_2(m_r, m_fe, m_K, m_p_new, m_x_new);
			}
			else if(m_system_1)
			{
				m_system_1(m_r, m_K, m_x_new);
				for(uint32_t i = 0; i < m_size; i++) m_r[i] = m_p_new * m_fe[i] - m_r[i];
			}
			else
			{
				m_residue(m_r, m_p_new, m_x_new);
				m_tangent_2(m_K, m_p_new, m_x_new);
				m_tangent_1(m_fe, m_p_new, m_x_new);
			}
		}
		void NewtonRaphson::predictor(void)
		{
			//predictor
			if(!solve(m_K, m_r, m_dxr) || !solve(m_K, m_fe, m_dxt))
			{
				if(!m_silent) printf("Unable to decompose stiffness Matrix in predictor!\n");
			}
			//continuation
			if(m_step != 1)
			{
				m_dp = m_continuation.predictor() / (1 << m_attempt);
			}
			for(uint32_t i = 0; i < m_size; i++) m_dx[i] = m_dxr[i] + m_dp * m_dxt[i];
		}
		void NewtonRaphson::corrector(void)
		{
			//corrector
			if(!solve(m_K, m_r, m_ddxr) || !solve(m_K, m_fe, m_ddxt))
			{
				if(!m_silent) printf("Unable to decompose stiffness Matrix in corrector!\n");
			}
			m_ddp = m_continuation.corrector();
			//update
			m_dp += m_ddp;
			for(uint32_t i = 0; i < m_size; i++) m_dx[i] += m_ddxr[i] + m_ddp * m_ddxt[i];
		}
	}
}