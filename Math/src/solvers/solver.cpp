//std
#include <cstdio>
#include <cstring>

//math
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/solver.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		solver::solver(void) : 
			m_silent(false), 
			m_equilibrium(false),
			m_size(1), m_watch_dof(0),
			m_step(0), m_attempt(0), m_iteration(0),
			m_step_max(100), m_attempt_max(5), m_iteration_max(10),
			m_K(nullptr), m_C(nullptr), m_M(nullptr),
			m_r(nullptr), m_fi(nullptr), m_fe(nullptr), 
			m_dxr(nullptr), m_dxt(nullptr), m_ddxr(nullptr), m_ddxt(nullptr),
			m_x_old(nullptr), m_x_new(nullptr), m_x_data(nullptr), m_dx(nullptr),
			m_v_old(nullptr), m_v_new(nullptr), m_v_data(nullptr), m_dv(nullptr),
			m_a_old(nullptr), m_a_new(nullptr), m_a_data(nullptr), m_da(nullptr),
			m_t_old(0), m_t_new(0), m_t_data(nullptr), m_dt(0), m_t_max(1.00e+00),
			m_p_old(0), m_p_new(0), m_p_data(nullptr), m_dp(0), m_dp0(1.00e-01), m_ddp(0)
		{
			return;
		}

		//destructor
		solver::~solver(void)
		{
			solver::cleanup();
		}

		//data
		void solver::save(const char* path) const
		{
			//file
			FILE* file = fopen(path, "w");
			const uint32_t ss = state_set();
			//write
			for(uint32_t i = 0; i < m_step; i++)
			{
				for(uint32_t j = 0; j < m_size; j++)
				{
					if(ss & uint32_t(state::x)) fprintf(file, "%+.6e ", m_x_data[j + m_size * i]);
					if(ss & uint32_t(state::v)) fprintf(file, "%+.6e ", m_v_data[j + m_size * i]);
					if(ss & uint32_t(state::a)) fprintf(file, "%+.6e ", m_a_data[j + m_size * i]);
				}
				if(ss & uint32_t(state::t)) fprintf(file, "%+.6e ", m_t_data[i]);
				if(ss & uint32_t(state::p)) fprintf(file, "%+.6e ", m_p_data[i]);
				fprintf(file, "\n");
			}
			//close
			fclose(file);
		}

		//solve
		bool solver::stop(void)
		{
			return m_stop_criteria.stop() || (m_stop && m_stop());
		}
		void solver::apply(void)
		{
			//data
			const uint32_t ss = state_set();
			//apply
			for(uint32_t i = 0; i < m_size; i++)
			{
				if(ss & uint32_t(state::x)) m_x_new[i] = m_x_old[i] + m_dx[i];
				if(ss & uint32_t(state::v)) m_v_new[i] = m_v_old[i] + m_dv[i];
				if(ss & uint32_t(state::a)) m_a_new[i] = m_a_old[i] + m_da[i];
			}
			if(ss & uint32_t(state::t)) m_t_new = m_t_old + m_dt;
			if(ss & uint32_t(state::p)) m_p_new = m_p_old + m_dp;
		}
		void solver::print(void)
		{
			//data
			const uint32_t ss = state_set();
			//print
			if(!m_silent)
			{
				printf("step: %04d ", m_step);
				if(ss & uint32_t(state::t)) printf("time: %+.6e ", m_t_new);
				if(ss & uint32_t(state::p)) printf("load: %+.6e ", m_p_new);
				if(ss & uint32_t(state::x)) printf("state: %+.6e ", m_x_new[m_watch_dof]);
				if(ss & uint32_t(state::v)) printf("velocity: %+.6e ", m_v_new[m_watch_dof]);
				if(ss & uint32_t(state::a)) printf("acceleration: %+.6e ", m_a_new[m_watch_dof]);
				printf("\n");
			}
			if(m_interface) m_interface(m_step);
		}
		void solver::setup(void)
		{
			//data
			compute();
			const uint32_t ss = state_set();
			//setup
			m_step = 0;
			m_dp = m_dp0;
			m_p_old = m_p_new;
			m_t_old = m_t_new;
			m_dt = m_t_max / m_step_max;
			m_convergence.m_solver = this;
			m_continuation.m_solver = this;
			m_stop_criteria.m_solver = this;
			if(ss & uint32_t(state::x)) memcpy(m_x_old, m_x_new, m_size * sizeof(double));
			if(ss & uint32_t(state::v)) memcpy(m_v_old, m_v_new, m_size * sizeof(double));
			if(ss & uint32_t(state::a)) memcpy(m_a_old, m_a_new, m_size * sizeof(double));
		}
		void solver::record(void)
		{
			//data
			const uint32_t ss = state_set();
			//record
			if(m_record) m_record();
			for(uint32_t i = 0; i < m_size; i++)
			{
				if(ss & uint32_t(state::x)) m_x_data[m_step * m_size + i] = m_x_new[i];
				if(ss & uint32_t(state::v)) m_v_data[m_step * m_size + i] = m_v_new[i];
				if(ss & uint32_t(state::a)) m_a_data[m_step * m_size + i] = m_a_new[i];
			}
			if(ss & uint32_t(state::t)) m_t_data[m_step] = m_t_new;
			if(ss & uint32_t(state::p)) m_p_data[m_step] = m_p_new;
		}
		void solver::update(void)
		{
			//data
			const uint32_t ss = state_set();
			//update
			if(m_update) m_update();
			if(ss & uint32_t(state::t)) m_t_old = m_t_new;
			if(ss & uint32_t(state::p)) m_p_old = m_p_new;
			if(ss & uint32_t(state::x)) memcpy(m_x_old, m_x_new, m_size * sizeof(double));
			if(ss & uint32_t(state::v)) memcpy(m_v_old, m_v_new, m_size * sizeof(double));
			if(ss & uint32_t(state::a)) memcpy(m_a_old, m_a_new, m_size * sizeof(double));
		}
		void solver::restore(void)
		{
			//data
			const uint32_t ss = state_set();
			//update
			if(m_update) m_update();
			if(ss & uint32_t(state::t)) m_t_new = m_t_old;
			if(ss & uint32_t(state::p)) m_p_new = m_p_old;
			if(ss & uint32_t(state::x)) memcpy(m_x_new, m_x_old, m_size * sizeof(double));
			if(ss & uint32_t(state::v)) memcpy(m_v_new, m_v_old, m_size * sizeof(double));
			if(ss & uint32_t(state::a)) memcpy(m_a_new, m_a_old, m_size * sizeof(double));
		}
		bool solver::equilibrium(void)
		{
			return m_equilibrium = m_convergence.check();
		}

		//allocate
		void solver::allocate_state(void)
		{
			const uint32_t ss = state_set();
			if(ss & uint32_t(state::x)) m_dx = new double[m_size];
			if(ss & uint32_t(state::v)) m_dv = new double[m_size];
			if(ss & uint32_t(state::a)) m_da = new double[m_size];
			if(ss & uint32_t(state::x)) m_dxr = new double[m_size];
			if(ss & uint32_t(state::x)) m_dxt = new double[m_size];
			if(ss & uint32_t(state::x)) m_ddxr = new double[m_size];
			if(ss & uint32_t(state::x)) m_ddxt = new double[m_size];
			if(ss & uint32_t(state::x)) m_x_old = new double[m_size];
			if(ss & uint32_t(state::x)) m_x_new = new double[m_size];
			if(ss & uint32_t(state::v)) m_v_old = new double[m_size];
			if(ss & uint32_t(state::v)) m_v_new = new double[m_size];
			if(ss & uint32_t(state::a)) m_a_old = new double[m_size];
			if(ss & uint32_t(state::a)) m_a_new = new double[m_size];
			if(ss & uint32_t(state::t)) m_t_data = new double[m_step_max + 1];
			if(ss & uint32_t(state::p)) m_p_data = new double[m_step_max + 1];
			if(ss & uint32_t(state::x)) m_x_data = new double[m_size * (m_step_max + 1)];
			if(ss & uint32_t(state::v)) m_v_data = new double[m_size * (m_step_max + 1)];
			if(ss & uint32_t(state::a)) m_a_data = new double[m_size * (m_step_max + 1)];
		}
		void solver::allocate_forces(void)
		{
			const uint32_t fs = force_set();
			if(fs & uint32_t(force::r)) m_r = new double[m_size];
			if(fs & uint32_t(force::fi)) m_fi = new double[m_size];
			if(fs & uint32_t(force::fe)) m_fe = new double[m_size];
		}
		void solver::allocate_tangents(void)
		{
			const uint32_t ts = tangent_set();
			if(ts & uint32_t(tangent::K)) m_K = new double[m_size * m_size];
			if(ts & uint32_t(tangent::C)) m_C = new double[m_size * m_size];
			if(ts & uint32_t(tangent::M)) m_M = new double[m_size * m_size];
		}

		//solve
		void solver::step(void)
		{
			for(m_attempt = 0; m_attempt < m_attempt_max; m_attempt++)
			{
				predictor();
				for(m_iteration = 0; m_iteration < m_iteration_max; m_iteration++)
				{
					apply();
					compute();
					if(equilibrium()) break; else corrector();
				}
				if(m_equilibrium) break;
				restore();
			}
			update();
			record();
		}
		void solver::solve(void)
		{
			check();
			setup();
			print();
			record();
			for(m_step = 1; !stop(); m_step++)
			{
				step();
				print();
				if(!m_equilibrium)
				{
					if(!m_silent) printf("Solver failed in step %d!\n", m_step);
					break;
				}
			}
		}
		void solver::cleanup(void)
		{
			//data
			double** data[] = {
				&m_K, &m_C, &m_M, 
				&m_r, &m_fi, &m_fe,
				&m_dxr, &m_dxt, &m_ddxr, &m_ddxt,
				&m_x_old, &m_x_new, &m_x_data, &m_dx,
				&m_v_old, &m_v_new, &m_v_data, &m_dv,
				&m_a_old, &m_a_new, &m_a_data, &m_da,
				&m_t_data, &m_p_data
			};
			//cleanup
			for(double** ptr : data)
			{
				delete[] *ptr;
				*ptr = nullptr;
			}
		}
		void solver::allocate(void)
		{
			allocate_state();
			allocate_forces();
			allocate_tangents();
		}
		void solver::allocate(uint32_t size)
		{
			m_size = size;
			allocate();
		}
	}
}