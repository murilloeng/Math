//std
#include <cstdio>
#include <cstring>

//Math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Linear/Sparse.hpp"
#include "Math/inc/Solvers/Solver.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		Solver::Solver(void) : 
			m_silent{false},
			m_rows_map{nullptr}, m_cols_map{nullptr},
			m_size{1}, m_watch_dof{0}, 
			m_K{nullptr}, m_C{nullptr}, m_M{nullptr},
			m_r{nullptr}, m_fi{nullptr}, m_fe{nullptr}, 
			m_dxr{nullptr}, m_dxt{nullptr}, m_ddxr{nullptr}, m_ddxt{nullptr},
			m_x_old{nullptr}, m_x_new{nullptr}, m_dx{nullptr},
			m_v_old{nullptr}, m_v_new{nullptr}, m_dv{nullptr},
			m_a_old{nullptr}, m_a_new{nullptr}, m_da{nullptr},
			m_p_old{0}, m_p_new{0}, m_dp{0}, m_dp0{1.00e-02}, m_ddp{0},
			m_t_old{0}, m_t_new{0}, m_dt{0}, m_t_min{0.00e+00}, m_t_max{1.00e+00}
		{
			return;
		}

		//destructor
		Solver::~Solver(void)
		{
			//data
			const double* data[] = {
				m_K, m_C, m_M, 
				m_r, m_fi, m_fe,
				m_dxr, m_dxt, m_ddxr, m_ddxt,
				m_x_old, m_x_new, m_dx, m_v_old, m_v_new, m_dv, m_a_old, m_a_new, m_da
			};
			//delete
			for(const double* ptr : data)
			{
				delete[] ptr;
			}
		}

		//serialization
		void Solver::save(const char*) const
		{
			return;
		}

		//solve
		bool Solver::stop(void)
		{
			return m_stop && m_stop();
		}
		void Solver::apply(void)
		{
			//data
			const uint32_t ss = state_set();
			//apply
			for(uint32_t i = 0; i < m_size; i++)
			{
				if(ss & uint32_t(State::x)) m_x_new[i] = m_x_old[i] + m_dx[i];
				if(ss & uint32_t(State::v)) m_v_new[i] = m_v_old[i] + m_dv[i];
				if(ss & uint32_t(State::a)) m_a_new[i] = m_a_old[i] + m_da[i];
			}
			if(ss & uint32_t(State::t)) m_t_new = m_t_old + m_dt;
			if(ss & uint32_t(State::p)) m_p_new = m_p_old + m_dp;
		}
		void Solver::print(void)
		{
			return;
		}
		void Solver::setup(void)
		{
			m_dp = m_dp0;
			m_p_old = m_p_new;
			m_t_old = m_t_new = m_t_min;
			if(state_set() & uint32_t(State::x)) memcpy(m_x_old, m_x_new, m_size * sizeof(double));
			if(state_set() & uint32_t(State::v)) memcpy(m_v_old, m_v_new, m_size * sizeof(double));
			if(state_set() & uint32_t(State::a)) memcpy(m_a_old, m_a_new, m_size * sizeof(double));
		}
		void Solver::record(void)
		{
			return;
		}
		void Solver::update(void)
		{
			//data
			const uint32_t ss = state_set();
			//update
			if(m_update) m_update();
			if(ss & uint32_t(State::t)) m_t_old = m_t_new;
			if(ss & uint32_t(State::p)) m_p_old = m_p_new;
			if(ss & uint32_t(State::x)) memcpy(m_x_old, m_x_new, m_size * sizeof(double));
			if(ss & uint32_t(State::v)) memcpy(m_v_old, m_v_new, m_size * sizeof(double));
			if(ss & uint32_t(State::a)) memcpy(m_a_old, m_a_new, m_size * sizeof(double));
		}
		void Solver::restore(void)
		{
			//data
			const uint32_t ss = state_set();
			//update
			if(m_restore) m_restore();
			if(ss & uint32_t(State::t)) m_t_new = m_t_old;
			if(ss & uint32_t(State::p)) m_p_new = m_p_old;
			if(ss & uint32_t(State::x)) memcpy(m_x_new, m_x_old, m_size * sizeof(double));
			if(ss & uint32_t(State::v)) memcpy(m_v_new, m_v_old, m_size * sizeof(double));
			if(ss & uint32_t(State::a)) memcpy(m_a_new, m_a_old, m_size * sizeof(double));
		}

		//allocate
		void Solver::allocate_state(void)
		{
			const uint32_t ss = state_set();
			if(ss & uint32_t(State::x)) m_dx = new double[m_size];
			if(ss & uint32_t(State::v)) m_dv = new double[m_size];
			if(ss & uint32_t(State::a)) m_da = new double[m_size];
			if(ss & uint32_t(State::x)) m_dxr = new double[m_size];
			if(ss & uint32_t(State::x)) m_dxt = new double[m_size];
			if(ss & uint32_t(State::x)) m_ddxr = new double[m_size];
			if(ss & uint32_t(State::x)) m_ddxt = new double[m_size];
			if(ss & uint32_t(State::x)) m_x_old = new double[m_size];
			if(ss & uint32_t(State::x)) m_x_new = new double[m_size];
			if(ss & uint32_t(State::v)) m_v_old = new double[m_size];
			if(ss & uint32_t(State::v)) m_v_new = new double[m_size];
			if(ss & uint32_t(State::a)) m_a_old = new double[m_size];
			if(ss & uint32_t(State::a)) m_a_new = new double[m_size];
		}
		void Solver::allocate_forces(void)
		{
			const uint32_t fs = force_set();
			if(fs & uint32_t(Force::r)) m_r = new double[m_size];
			if(fs & uint32_t(Force::fi)) m_fi = new double[m_size];
			if(fs & uint32_t(Force::fe)) m_fe = new double[m_size];
		}
		void Solver::allocate_tangents(void)
		{
			const uint32_t ts = tangent_set();
			if(ts & uint32_t(Tangent::K)) m_K = new double[m_size * m_size];
			if(ts & uint32_t(Tangent::C)) m_C = new double[m_size * m_size];
			if(ts & uint32_t(Tangent::M)) m_M = new double[m_size * m_size];
		}

		//solve
		bool Solver::solve(const double* K, const double* f, double* x) const
		{
			if(m_rows_map == nullptr || m_cols_map == nullptr)
			{
				//data
				math::Vector xm(x, m_size);
				const math::Vector fm(f, m_size);
				const math::Matrix Km(K, m_size, m_size);
				//solve
				return Km.solve(xm, fm);
			}
			else
			{
				//data
				math::Vector xm(x, m_size);
				const math::Vector fm(f, m_size);
				const math::Sparse Km(K, m_rows_map, m_cols_map, m_size, m_size);
				//solve
				return Km.solve(xm, fm);
			}
		}

		//solve
		void Solver::step(void)
		{
			// for(m_attempt = 0; m_attempt < m_attempt_max; m_attempt++)
			// {
			// 	predictor();
			// 	for(m_iteration = 0; m_iteration < m_iteration_max; m_iteration++)
			// 	{
			// 		apply();
			// 		compute();
			// 		if(equilibrium()) break; else corrector();
			// 	}
			// 	if(m_equilibrium) break;
			// 	restore();
			// }
			// update();
			// record();
		}
		void Solver::solve(void)
		{
			// check();
			// setup();
			// print();
			// record();
			// compute();
			// for(m_step = 1; !stop(); m_step++)
			// {
			// 	step();
			// 	print();
			// 	if(!m_equilibrium)
			// 	{
			// 		if(!m_silent) printf("Solver failed in step %d!\n", m_step);
			// 		break;
			// 	}
			// }
		}
		void Solver::cleanup(void)
		{
			//data
			double** data[] = {
				&m_K, &m_C, &m_M, 
				&m_r, &m_fi, &m_fe,
				&m_dxr, &m_dxt, &m_ddxr, &m_ddxt,
				&m_x_old, &m_x_new, &m_dx, &m_v_old, &m_v_new, &m_dv, &m_a_old, &m_a_new, &m_da
			};
			//cleanup
			for(double** ptr : data)
			{
				delete[] *ptr;
				*ptr = nullptr;
			}
		}
		void Solver::allocate(void)
		{
			allocate_state();
			allocate_forces();
			allocate_tangents();
		}
		void Solver::allocate(uint32_t size)
		{
			m_size = size;
			allocate_state();
			allocate_forces();
			allocate_tangents();
		}
	}
}