//std
#include <cmath>
#include <cstdio>
#include <cstring>

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/linear/matrix.hpp"
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/newton_raphson.hpp"

namespace math
{
	namespace solvers
	{
		//constructors
		newton_raphson::newton_raphson(void) : 
			m_silent(false), 
			m_watch_dof(0), m_size(1), m_step_max(100), m_attempt_max(5), m_iteration_max(10),
			m_dp0(0.01), m_tolerance(1e-5), m_p_old(0), m_p_new(0), m_p_data(nullptr),
			m_x_old(nullptr), m_x_new(nullptr), m_x_data(nullptr),
			m_r(nullptr), m_g(nullptr), m_K(nullptr), 
			m_dx(nullptr), m_dxr(nullptr), m_dxt(nullptr), m_ddxr(nullptr), m_ddxt(nullptr)
		{
			return;
		}

		//destructor
		newton_raphson::~newton_raphson(void)
		{
			cleanup();
		}

		//data
		void newton_raphson::size(uint32_t size)
		{
			//data
			m_size = size;
			//memory
			cleanup();
			allocate();
		}
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

		//solve
		bool newton_raphson::stop(void)
		{
			return (m_stop && m_stop()) || m_stop_criteria.stop();
		}
		bool newton_raphson::check(void)
		{
			return m_system_2 || m_system_1 || (m_residue && m_tangent_1 && m_tangent_2);
		}
		void newton_raphson::apply(void)
		{
			m_p_new = m_p_old + m_dp;
			for(uint32_t i = 0; i < m_size; i++)
			{
				m_x_new[i] = m_x_old[i] + m_dx[i];
			}
			compute();
		}
		void newton_raphson::print(void)
		{
			if(!m_silent)
			{
				printf("step: %04d load: %+.6e state: %+.6e\n", m_step, m_p_new, m_x_new[m_watch_dof]);
			}
			if(m_interface) m_interface(m_step);
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
		void newton_raphson::update(void)
		{
			m_p_old = m_p_new;
			if(m_update) m_update();
			memcpy(m_x_old, m_x_new, m_size * sizeof(double));
		}
		void newton_raphson::record(void)
		{
			if(m_record) m_record();
			m_p_data[m_step] = m_p_new;
			for(uint32_t i = 0; i < m_size; i++)
			{
				m_x_data[m_step * m_size + i] = m_x_new[i];
			}
		}
		void newton_raphson::restore(void)
		{
			m_p_new = m_p_old;
			if(m_restore) m_restore();
			memcpy(m_x_new, m_x_old, m_size * sizeof(double));
		}
		void newton_raphson::compute(void)
		{
			if(m_system_2)
			{
				m_system_2(m_r, m_g, m_K, m_p_new, m_x_new);
			}
			else if(m_system_1)
			{
				m_system_1(m_r, m_K, m_x_new);
				for(uint32_t i = 0; i < m_size; i++) m_r[i] = m_p_new * m_g[i] - m_r[i];
			}
			else
			{
				m_residue(m_r, m_p_new, m_x_new);
				m_tangent_1(m_g, m_p_new, m_x_new);
				m_tangent_2(m_K, m_p_new, m_x_new);
			}
		}
		bool newton_raphson::equilibrium(void)
		{
			const vector r(m_r, m_size), g(m_g, m_size);
			return m_equilibrium = r.norm() < m_tolerance * g.norm();
		}
		void newton_raphson::predictor(void)
		{
			//data
			const matrix K(m_K, m_size, m_size);
			const vector r(m_r, m_size), g(m_g, m_size);
			vector dx(m_dx, m_size), dxr(m_dxr, m_size), dxt(m_dxt, m_size);
			//predictor
			if(!K.solve(dxr, r) || !K.solve(dxt, g))
			{
				printf("Unable to decompose stiffness matrix in predictor!\n");
			}
			load_predictor();
			for(uint32_t i = 0; i < m_size; i++) dx[i] = dxr[i] + m_dp * dxt[i];
			//apply
			apply();
		}
		void newton_raphson::corrector(void)
		{
			const matrix K(m_K, m_size, m_size);
			const vector r(m_r, m_size), g(m_g, m_size);
			vector ddxr(m_ddxr, m_size), ddxt(m_ddxt, m_size);
			for(m_iteration = 0; m_iteration < m_iteration_max; m_iteration++)
			{
				//check
				if(equilibrium()) break;
				//corrector
				if(!K.solve(ddxr, r) || !K.solve(ddxt, g))
				{
					printf("Unable to decompose stiffness matrix in corrector!\n");
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

		//solve
		void newton_raphson::step(void)
		{
			for(m_attempt = 0; m_attempt < m_attempt_max; m_attempt++)
			{
				predictor();
				corrector();
				if(m_equilibrium) break;
				restore();
			}
			update();
			record();
		}
		void newton_raphson::solve(void)
		{
			if(check())
			{
				setup();
				print();
				record();
				for(m_step = 1; !stop(); m_step++)
				{
					step();
					print();
					if(!m_equilibrium) break;
				}
			}
			else
			{
				printf("Error: Newton-Raphson solver failed because system methods are not set!\n");
			}
		}
		void newton_raphson::cleanup(void)
		{
			//data
			double** data[] = {
				&m_p_data, &m_x_old, &m_x_new, &m_x_data, 
				&m_r, &m_g, &m_K, &m_dx, &m_dxr, &m_dxt, &m_ddxr, &m_ddxt
			};
			//delete
			for(double** ptr : data)
			{
				delete[] *ptr;
				*ptr = nullptr;
			}
		}
		void newton_raphson::allocate(void)
		{
			m_r = new double[m_size];
			m_g = new double[m_size];
			m_dx = new double[m_size];
			m_dxr = new double[m_size];
			m_dxt = new double[m_size];
			m_ddxr = new double[m_size];
			m_ddxt = new double[m_size];
			m_x_old = new double[m_size];
			m_x_new = new double[m_size];
			m_K = new double[m_size * m_size];
			m_p_data = new double[m_step_max + 1];
			m_x_data = new double[m_size * (m_step_max + 1)];
		}
	}
}