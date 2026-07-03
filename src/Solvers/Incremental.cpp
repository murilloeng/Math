//std
#include <cstdio>

//Math
#include "Math/inc/Solvers/Implicit.hpp"
#include "Math/inc/Solvers/Incremental.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		Incremental::Incremental(void) :
			m_step{0}, m_step_max{100}, m_stop_criteria{this},
			m_x_data{nullptr}, m_v_data{nullptr}, m_a_data{nullptr}, m_p_data{nullptr}, m_t_data{nullptr}
		{
			return;
		}

		//destructor
		Incremental::~Incremental(void)
		{
			//data
			double** data[] = {
				&m_x_data, &m_v_data, &m_a_data, &m_p_data, &m_t_data
			};
			//cleanup
			for(double** ptr : data)
			{
				delete[] *ptr;
				*ptr = nullptr;
			}
		}

		//solve
		void Incremental::step(void)
		{
			return;
		}
		void Incremental::solve(void)
		{
			check();
			setup();
			print();
			record();
			compute();
			for(m_step = 1; !stop(); m_step++)
			{
				step();
				print();
				if(!m_status)
				{
					if(!m_silent) printf("Solver failed in step %d!\n", m_step);
					break;
				}
			}
		}

		//serialization
		void Incremental::save(const char* path) const
		{
			FILE* file = fopen(path, "w");
			const uint32_t ss = state_set();
			for(uint32_t i = 0; i < m_step; i++)
			{
				for(uint32_t j = 0; j < m_size; j++)
				{
					if(ss & 1 << uint32_t(State::x)) fprintf(file, "%+.6e ", m_x_data[j + m_size * i]);
					if(ss & 1 << uint32_t(State::v)) fprintf(file, "%+.6e ", m_v_data[j + m_size * i]);
					if(ss & 1 << uint32_t(State::a)) fprintf(file, "%+.6e ", m_a_data[j + m_size * i]);
				}
				if(ss & 1 << uint32_t(State::p)) fprintf(file, "%+.6e ", m_p_data[i]);
				if(ss & 1 << uint32_t(State::t)) fprintf(file, "%+.6e ", m_t_data[i]);
				fprintf(file, "\n");
			}
			fclose(file);
		}

		//solve
		bool Incremental::stop(void)
		{
			return Solver::stop() && m_stop_criteria.stop();
		}
		void Incremental::print(void)
		{
			if(m_silent) return;
			printf("Step: %4d ", m_step);
		}
		void Incremental::setup(void)
		{
			m_step = 0;
			m_dt = (m_t_max - m_t_min) / m_step_max;
		}
		void Incremental::record(void)
		{
			//data
			const uint32_t ss = state_set();
			//callback
			if(m_callback_record) m_callback_record();
			//record
			for(uint32_t i = 0; i < m_size; i++)
			{
				if(ss & uint32_t(State::x)) m_x_data[m_step * m_size + i] = m_x_new[i];
				if(ss & uint32_t(State::v)) m_v_data[m_step * m_size + i] = m_v_new[i];
				if(ss & uint32_t(State::a)) m_a_data[m_step * m_size + i] = m_a_new[i];
			}
			if(ss & uint32_t(State::t)) m_t_data[m_step] = m_t_new;
			if(ss & uint32_t(State::p)) m_p_data[m_step] = m_p_new;
		}

		//allocate
		void Incremental::cleanup(void)
		{
			//data
			double** data[] = {
				&m_x_data, &m_v_data, &m_a_data, &m_p_data, &m_t_data
			};
			//cleanup
			Solver::cleanup();
			for(double** ptr : data)
			{
				delete[] *ptr;
				*ptr = nullptr;
			}
		}
		void Incremental::allocate_state(void)
		{
			Solver::allocate_state();
			const uint32_t ss = state_set();
			if(ss & uint32_t(State::t)) m_t_data = new double[m_step_max + 1];
			if(ss & uint32_t(State::p)) m_p_data = new double[m_step_max + 1];
			if(ss & uint32_t(State::x)) m_x_data = new double[m_size * (m_step_max + 1)];
			if(ss & uint32_t(State::v)) m_v_data = new double[m_size * (m_step_max + 1)];
			if(ss & uint32_t(State::a)) m_a_data = new double[m_size * (m_step_max + 1)];
		}
	}
}