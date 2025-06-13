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
		void newton_raphson::save(const char* path) const
		{
			FILE* file = fopen(path, "w");
			fprintf(file, "%04d\n", m_step);
			for(uint32_t i = 0; i < m_step; i++)
			{
				fprintf(file, "%+.6e ", m_p_data[i]);
				for(uint32_t j = 0; j < m_size; j++)
				{
					fprintf(file, "%+.6e ", m_x_data[j + m_size * i]);
				}
				fprintf(file, "\n");
			}
		}
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
				printf("Error: Newton-Raphson solver called with methods not set!\n");
				exit(EXIT_FAILURE);
			}
		}
		void newton_raphson::setup(void)
		{
			//data
			compute();
			m_step = 0;
			m_dp = m_dp0;
			m_p_old = m_p_new;
			memcpy(m_x_old, m_x_new, m_size * sizeof(double));
			//stop
			m_stop_criteria.m_solver = this;
			//convergence
			m_convergence.m_solver = this;
			//continuation
			m_continuation.m_dx = m_dx;
			m_continuation.m_dp = &m_dp;
			m_continuation.m_dxr = m_dxr;
			m_continuation.m_dxt = m_dxt;
			m_continuation.m_ddxr = m_ddxr;
			m_continuation.m_ddxt = m_ddxt;
			m_continuation.m_size = m_size;
			m_continuation.m_index = m_watch_dof;
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
			load_predictor();
			for(uint32_t i = 0; i < m_size; i++) dx[i] = dxr[i] + m_dp * dxt[i];
			//apply
			apply();
		}
		void newton_raphson::corrector(void)
		{
			const matrix K(m_K, m_size, m_size);
			const vector r(m_r, m_size), fe(m_fe, m_size);
			vector ddxr(m_ddxr, m_size), ddxt(m_ddxt, m_size);
			for(m_iteration = 0; m_iteration < m_iteration_max; m_iteration++)
			{
				//check
				if(equilibrium()) break;
				//corrector
				if(!K.solve(ddxr, r) || !K.solve(ddxt, fe))
				{
					if(!m_silent) printf("Unable to decompose stiffness matrix in corrector!\n");
				}
				load_corrector();
				//update
				m_dp += m_ddp;
				for(uint32_t i = 0; i < m_size; i++) m_dx[i] += m_ddxr[i] + m_ddp * m_ddxt[i];
				//apply
				apply();
			}
		}
		void newton_raphson::load_predictor(void)
		{
			if(m_step != 1)
			{
				m_dp = m_continuation.predictor() / (1 << m_attempt);
			}
		}
		void newton_raphson::load_corrector(void)
		{
			m_ddp = m_continuation.corrector();
		}
	}
}