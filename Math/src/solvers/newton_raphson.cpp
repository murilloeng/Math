//std
#include <cmath>
#include <cstdio>
#include <cstring>

//math
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/newton_raphson.hpp"

namespace math
{
	namespace solvers
	{
		//constructors
		newton_raphson::newton_raphson(void)
		{
			return;
		}

		//destructor
		newton_raphson::~newton_raphson(void)
		{
			cleanup();
		}

		//data
		uint32_t newton_raphson::state_set(void) const
		{
			return uint32_t(state::x) | uint32_t(state::p);
		}
		uint32_t newton_raphson::force_set(void) const
		{
			return uint32_t(force::r) | uint32_t(force::fi) | uint32_t(force::fe);
		}
		uint32_t newton_raphson::tangent_set(void) const
		{
			return uint32_t(tangent::K);
		}

		//solve
		void newton_raphson::check(void)
		{
			if(!m_system_1 && !m_system_2 && !(m_residue && m_tangent_1 && m_tangent_2))
			{
				printf("Error: Newton-Raphson solver called with at least one method not set!\n");
				exit(EXIT_FAILURE);
			}
		}
		void newton_raphson::compute(void)
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
		void newton_raphson::predictor(void)
		{
			//data
			const matrix K(m_K, m_size, m_size);
			const vector r(m_r, m_size), fe(m_fe, m_size);
			vector dx(m_dx, m_size), dxr(m_dxr, m_size), dxt(m_dxt, m_size);
			//predictor
			if(!K.solve(dxr, r) || !K.solve(dxt, fe))
			{
				if(!m_silent) printf("Unable to decompose stiffness matrix in predictor!\n");
			}
			//continuation
			if(m_step != 1)
			{
				m_dp = m_continuation.predictor() / (1 << m_attempt);
			}
			for(uint32_t i = 0; i < m_size; i++) dx[i] = dxr[i] + m_dp * dxt[i];
		}
		void newton_raphson::corrector(void)
		{
			//data
			const matrix K(m_K, m_size, m_size);
			const vector r(m_r, m_size), fe(m_fe, m_size);
			vector ddxr(m_ddxr, m_size), ddxt(m_ddxt, m_size);
			//corrector
			if(!K.solve(ddxr, r) || !K.solve(ddxt, fe))
			{
				if(!m_silent) printf("Unable to decompose stiffness matrix in corrector!\n");
			}
			//continuation
			m_ddp = m_continuation.corrector();
			//update
			m_dp += m_ddp;
			for(uint32_t i = 0; i < m_size; i++) m_dx[i] += m_ddxr[i] + m_ddp * m_ddxt[i];
		}
	}
}