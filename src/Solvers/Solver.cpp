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
			m_silent{true},
			m_status{false},
			m_rows_map{nullptr},
			m_cols_map{nullptr},
			m_size{1}, m_watch_dof{0},
			m_K{nullptr}, m_C{nullptr}, m_M{nullptr},
			m_r{nullptr}, m_fi{nullptr}, m_fe{nullptr},
			m_x_old{nullptr}, m_x_new{nullptr}, m_dx{nullptr},
			m_v_old{nullptr}, m_v_new{nullptr}, m_dv{nullptr},
			m_a_old{nullptr}, m_a_new{nullptr}, m_da{nullptr},
			m_dxr{nullptr}, m_dxt{nullptr}, m_ddxr{nullptr}, m_ddxt{nullptr},
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
				m_x_old, m_x_new, m_dx,
				m_v_old, m_v_new, m_dv,
				m_a_old, m_a_new, m_da,
				m_dxr, m_dxt, m_ddxr, m_ddxt
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

		//data
		bool Solver::silent(void) const
		{
			return m_silent;
		}
		bool Solver::silent(bool silent)
		{
			return m_silent = silent;
		}

		bool Solver::status(void) const
		{
			return m_status;
		}

		uint32_t Solver::size(void) const
		{
			return m_size;
		}
		uint32_t Solver::size(uint32_t size)
		{
			return m_size = size;
		}

		uint32_t Solver::watch_dof(void) const
		{
			return m_watch_dof;
		}
		uint32_t Solver::watch_dof(uint32_t watch_dof)
		{
			return m_watch_dof = watch_dof;
		}

		int32_t* Solver::rows_map(void) const
		{
			return m_rows_map;
		}
		int32_t* Solver::cols_map(void) const
		{
			return m_cols_map;
		}
		int32_t* Solver::rows_map(int32_t* rows_map)
		{
			return m_rows_map = rows_map;
		}
		int32_t* Solver::cols_map(int32_t* cols_map)
		{
			return m_cols_map = cols_map;
		}

		double Solver::time_min(void) const
		{
			return m_t_min;
		}
		double Solver::time_max(void) const
		{
			return m_t_max;
		}
		double Solver::time_min(double t_min)
		{
			return m_t_min = t_min;
		}
		double Solver::time_max(double t_max)
		{
			return m_t_max = t_max;
		}
		double Solver::time_old(void) const
		{
			return m_t_old;
		}
		double Solver::time_new(void) const
		{
			return m_t_new;
		}

		double Solver::step_size(void) const
		{
			return m_dp0;
		}
		double Solver::step_size(double dp0)
		{
			return m_dp0 = dp0;
		}

		double* Solver::state_old(void)
		{
			return m_x_old;
		}
		double* Solver::state_new(void)
		{
			return m_x_new;
		}
		double Solver::state_old(uint32_t index)
		{
			return m_x_old[index];
		}
		double Solver::state_new(uint32_t index)
		{
			return m_x_new[index];
		}
		double Solver::state_old(uint32_t index, double x)
		{
			return m_x_old[index] = x;
		}
		double Solver::state_new(uint32_t index, double x)
		{
			return m_x_new[index] = x;
		}

		double* Solver::velocity_old(void)
		{
			return m_v_old;
		}
		double* Solver::velocity_new(void)
		{
			return m_v_new;
		}
		double Solver::velocity_old(uint32_t index)
		{
			return m_v_old[index];
		}
		double Solver::velocity_new(uint32_t index)
		{
			return m_v_new[index];
		}
		double Solver::velocity_old(uint32_t index, double v)
		{
			return m_v_old[index] = v;
		}
		double Solver::velocity_new(uint32_t index, double v)
		{
			return m_v_new[index] = v;
		}

		double* Solver::acceleration_old(void)
		{
			return m_a_old;
		}
		double* Solver::acceleration_new(void)
		{
			return m_a_new;
		}
		double Solver::acceleration_old(uint32_t index)
		{
			return m_a_old[index];
		}
		double Solver::acceleration_new(uint32_t index)
		{
			return m_a_new[index];
		}
		double Solver::acceleration_old(uint32_t index, double a)
		{
			return m_a_old[index] = a;
		}
		double Solver::acceleration_new(uint32_t index, double a)
		{
			return m_a_new[index] = a;
		}

		//data
		void Solver::cleanup(void)
		{
			//data
			double** data[] = {
				&m_K, &m_C, &m_M,
				&m_r, &m_fi, &m_fe,
				&m_x_old, &m_x_new, &m_dx,
				&m_v_old, &m_v_new, &m_dv,
				&m_a_old, &m_a_new, &m_da,
				&m_dxr, &m_dxt, &m_ddxr, &m_ddxt
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

		//solve
		bool Solver::stop(void)
		{
			return m_callback_stop && m_callback_stop();
		}
		void Solver::check(void)
		{
			return;
		}
		void Solver::apply(void)
		{
			//data
			const uint32_t ss = state_set();
			//apply
			for(uint32_t i = 0; i < m_size; i++)
			{
				if(ss & 1 << uint32_t(State::x)) m_x_new[i] = m_x_old[i] + m_dx[i];
				if(ss & 1 << uint32_t(State::v)) m_v_new[i] = m_v_old[i] + m_dv[i];
				if(ss & 1 << uint32_t(State::a)) m_a_new[i] = m_a_old[i] + m_da[i];
			}
			if(ss & 1 << uint32_t(State::t)) m_t_new = m_t_old + m_dt;
			if(ss & 1 << uint32_t(State::p)) m_p_new = m_p_old + m_dp;
		}
		void Solver::print(void)
		{
			//data
			if(m_silent) return;
			const uint32_t ss = state_set();
			//print
			if(ss & 1 << uint32_t(State::t)) printf("Time: %+.6e ", m_t_new);
			if(ss & 1 << uint32_t(State::p)) printf("Load: %+.6e ", m_p_new);
			if(ss & 1 << uint32_t(State::x)) printf("State: %+.6e ", m_x_new[m_watch_dof]);
			if(ss & 1 << uint32_t(State::v)) printf("Velocity: %+.6e ", m_v_new[m_watch_dof]);
			if(ss & 1 << uint32_t(State::a)) printf("Acceleration: %+.6e ", m_a_new[m_watch_dof]);
			printf("\n");
		}
		void Solver::setup(void)
		{
			m_dp = m_dp0;
			m_status = false;
			m_p_new = m_p_old;
			m_t_new = m_t_old = m_t_min;
			const uint32_t ss = state_set();
			if(ss & 1 << uint32_t(State::x)) memcpy(m_x_new, m_x_old, m_size * sizeof(double));
			if(ss & 1 << uint32_t(State::v)) memcpy(m_v_new, m_v_old, m_size * sizeof(double));
			if(ss & 1 << uint32_t(State::a)) memcpy(m_a_new, m_a_old, m_size * sizeof(double));
		}
		void Solver::record(void)
		{
			if(m_callback_record) m_callback_record();
		}
		void Solver::update(void)
		{
			//data
			const uint32_t ss = state_set();
			//update
			if(m_callback_update) m_callback_update();
			if(ss & 1 << uint32_t(State::t)) m_t_old = m_t_new;
			if(ss & 1 << uint32_t(State::p)) m_p_old = m_p_new;
			if(ss & 1 << uint32_t(State::x)) memcpy(m_x_old, m_x_new, m_size * sizeof(double));
			if(ss & 1 << uint32_t(State::v)) memcpy(m_v_old, m_v_new, m_size * sizeof(double));
			if(ss & 1 << uint32_t(State::a)) memcpy(m_a_old, m_a_new, m_size * sizeof(double));
		}
		void Solver::restore(void)
		{
			//data
			const uint32_t ss = state_set();
			//update
			if(m_callback_restore) m_callback_restore();
			if(ss & 1 << uint32_t(State::t)) m_t_new = m_t_old;
			if(ss & 1 << uint32_t(State::p)) m_p_new = m_p_old;
			if(ss & 1 << uint32_t(State::x)) memcpy(m_x_new, m_x_old, m_size * sizeof(double));
			if(ss & 1 << uint32_t(State::v)) memcpy(m_v_new, m_v_old, m_size * sizeof(double));
			if(ss & 1 << uint32_t(State::a)) memcpy(m_a_new, m_a_old, m_size * sizeof(double));
		}
		void Solver::compute(void)
		{
			return;
		}
		void Solver::predictor(void)
		{
			return;
		}
		void Solver::corrector(void)
		{
			return;
		}

		//allocate
		void Solver::allocate_state(void)
		{
			const uint32_t ss = state_set();
			if(ss & 1 << uint32_t(State::x)) m_dx = new double[m_size];
			if(ss & 1 << uint32_t(State::v)) m_dv = new double[m_size];
			if(ss & 1 << uint32_t(State::a)) m_da = new double[m_size];
			if(ss & 1 << uint32_t(State::x)) m_dxr = new double[m_size];
			if(ss & 1 << uint32_t(State::x)) m_dxt = new double[m_size];
			if(ss & 1 << uint32_t(State::x)) m_ddxr = new double[m_size];
			if(ss & 1 << uint32_t(State::x)) m_ddxt = new double[m_size];
			if(ss & 1 << uint32_t(State::x)) m_x_old = new double[m_size];
			if(ss & 1 << uint32_t(State::x)) m_x_new = new double[m_size];
			if(ss & 1 << uint32_t(State::v)) m_v_old = new double[m_size];
			if(ss & 1 << uint32_t(State::v)) m_v_new = new double[m_size];
			if(ss & 1 << uint32_t(State::a)) m_a_old = new double[m_size];
			if(ss & 1 << uint32_t(State::a)) m_a_new = new double[m_size];
		}
		void Solver::allocate_forces(void)
		{
			const uint32_t fs = force_set();
			if(fs & 1 << uint32_t(Force::r)) m_r = new double[m_size];
			if(fs & 1 << uint32_t(Force::fi)) m_fi = new double[m_size];
			if(fs & 1 << uint32_t(Force::fe)) m_fe = new double[m_size];
		}
		void Solver::allocate_tangents(void)
		{
			const uint32_t ts = tangent_set();
			if(ts & 1 << uint32_t(Tangent::K)) m_K = new double[m_cols_map ? m_cols_map[m_size] : m_size * m_size];
			if(ts & 1 << uint32_t(Tangent::C)) m_C = new double[m_cols_map ? m_cols_map[m_size] : m_size * m_size];
			if(ts & 1 << uint32_t(Tangent::M)) m_M = new double[m_cols_map ? m_cols_map[m_size] : m_size * m_size];
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
	}
}