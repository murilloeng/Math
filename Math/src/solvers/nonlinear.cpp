//std
#include <cstdio>

//math
#include "Math/Math/inc/solvers/nonlinear.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		nonlinear::nonlinear(void) : 
			m_silent(true), 
			m_equilibrium(false),
			m_size(0), m_watch_dof(0),
			m_step(0), m_attempt(0), m_iteration(0),
			m_step_max(0), m_attempt_max(0), m_iteration_max(0),
			m_K(nullptr), m_C(nullptr), m_M(nullptr),
			m_r(nullptr), m_fi(nullptr), m_fe(nullptr), 
			m_dxr(nullptr), m_dxt(nullptr), m_ddxr(nullptr), m_ddxt(nullptr),
			m_x_old(nullptr), m_x_new(nullptr), m_x_data(nullptr), m_dx(nullptr),
			m_v_old(nullptr), m_v_new(nullptr), m_v_data(nullptr), m_dv(nullptr),
			m_a_old(nullptr), m_a_new(nullptr), m_a_data(nullptr), m_da(nullptr),
			m_p_old(0), m_p_new(0), m_p_data(nullptr), m_dp(0), m_dp0(0)
		{
			return;
		}

		//destructor
		nonlinear::~nonlinear(void)
		{
			nonlinear::cleanup();
		}

		//solve
		bool nonlinear::stop(void)
		{
			return (m_stop && m_stop()) || m_stop_criteria.stop();
			return true;
		}
		void nonlinear::check(void)
		{
			return;
		}
		void nonlinear::apply(void)
		{
			//data
			const uint32_t ss = state_set();
			//apply
			if(ss & uint32_t(state::p)) m_p_new = m_p_old + m_dp;
			for(uint32_t i = 0; i < m_size; i++)
			{
				if(ss & uint32_t(state::x)) m_x_new[i] = m_x_old[i] + m_dx[i];
				if(ss & uint32_t(state::v)) m_v_new[i] = m_v_old[i] + m_dv[i];
				if(ss & uint32_t(state::a)) m_a_new[i] = m_a_old[i] + m_da[i];
			}
			compute();
		}
		void nonlinear::print(void)
		{
			if(!m_silent)
			{
				printf("step: %04d state: %+.6e\n", m_step, m_x_new[m_watch_dof]);
			}
			if(m_interface) m_interface(m_step);
		}
		void nonlinear::setup(void)
		{
			// //data
			// compute();
			// m_step = 0;
			// m_dp = m_dp0;
			// m_p_old = m_p_new;
			// memcpy(m_x_old, m_x_new, m_size * sizeof(double));
			// //stop
			// m_stop_criteria.m_solver = this;
			// //convergence
			// m_convergence.m_r = m_r;
			// m_convergence.m_g = m_g;
			// //continuation
			// m_continuation.m_dx = m_dx;
			// m_continuation.m_dp = &m_dp;
			// m_continuation.m_dxr = m_dxr;
			// m_continuation.m_dxt = m_dxt;
			// m_continuation.m_ddxr = m_ddxr;
			// m_continuation.m_ddxt = m_ddxt;
			// m_continuation.m_size = m_size;
			// m_continuation.m_index = m_watch_dof;
		}
		void nonlinear::record(void)
		{
			return;
		}
		void nonlinear::update(void)
		{
			return;
		}
		void nonlinear::restore(void)
		{
			return;
		}
		void nonlinear::compute(void)
		{
			return;
		}
		void nonlinear::predictor(void)
		{
			return;
		}
		void nonlinear::corrector(void)
		{
			return;
		}
		bool nonlinear::equilibrium(void)
		{
			return true;
		}

		//solve
		void nonlinear::step(void)
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
		void nonlinear::solve(void)
		{
			check();
			setup();
			print();
			record();
			for(m_step = 0; !stop(); m_step++)
			{
				step();
				print();
				if(!m_equilibrium)
				{
					printf("Solver failed in step %d!\n", m_step);
					break;
				}
			}
		}
		void nonlinear::cleanup(void)
		{
			//data
			double** data[] = {
				&m_K, &m_C, &m_M, 
				&m_r, &m_fi, &m_fe, &m_p_data,
				&m_dxr, &m_dxt, &m_ddxr, &m_ddxt,
				&m_x_old, &m_x_new, &m_x_data, &m_dx,
				&m_v_old, &m_v_new, &m_v_data, &m_dv,
				&m_a_old, &m_a_new, &m_a_data, &m_da
			};
			//cleanup
			for(double** ptr : data)
			{
				delete[] *ptr;
				*ptr = nullptr;
			}
		}
		void nonlinear::allocate(void)
		{
			//data
			const uint32_t ss = state_set();
			const uint32_t fs = force_set();
			const uint32_t ts = tangent_set();
			//forces
			if(fs & uint32_t(force::r)) m_r = new double[m_size];
			if(fs & uint32_t(force::fi)) m_fi = new double[m_size];
			if(fs & uint32_t(force::fe)) m_fe = new double[m_size];
			//tangents
			if(ts & uint32_t(tangent::K)) m_K = new double[m_size * m_size];
			if(ts & uint32_t(tangent::C)) m_C = new double[m_size * m_size];
			if(ts & uint32_t(tangent::M)) m_M = new double[m_size * m_size];
			//state
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
			if(ss & uint32_t(state::p)) m_p_data = new double[m_step_max + 1];
			if(ss & uint32_t(state::x)) m_x_data = new double[m_size * (m_step_max + 1)];
			if(ss & uint32_t(state::v)) m_v_data = new double[m_size * (m_step_max + 1)];
			if(ss & uint32_t(state::a)) m_a_data = new double[m_size * (m_step_max + 1)];
		}
	}
}