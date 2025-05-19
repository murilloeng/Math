//std
#include <cfloat>

//Math
#include "Math/Math/inc/solvers/stop_criteria.hpp"
#include "Math/Math/inc/solvers/newton_raphson.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		stop_criteria::stop_criteria(void) : 
			m_stop(type::last), m_types(uint32_t(type::step_maximum)), 
			m_p_min(-DBL_MAX), m_p_max(+DBL_MAX), m_x_min(-DBL_MAX), m_x_max(+DBL_MAX)
		{
			return;
		}

		//destructor
		stop_criteria::~stop_criteria(void)
		{
			return;
		}

		//stop
		bool stop_criteria::stop(void)
		{
			//data
			bool(stop_criteria::*fun[])(void) const = {
				&stop_criteria::stop_step_maximum,
				&stop_criteria::stop_load_limit_minimum, &stop_criteria::stop_load_limit_maximum,
				&stop_criteria::stop_load_local_minimum, &stop_criteria::stop_load_local_maximum,
				&stop_criteria::stop_load_value_negative, &stop_criteria::stop_load_value_positive,
				&stop_criteria::stop_state_limit_minimum, &stop_criteria::stop_state_limit_maximum,
				&stop_criteria::stop_state_local_minimum, &stop_criteria::stop_state_local_maximum,
				&stop_criteria::stop_state_value_negative, &stop_criteria::stop_state_value_positive
			};
			//stop
			m_stop = type::last;
			for(uint32_t i = 0; 1U << i < uint32_t(type::last); i++)
			{
				if((i == 0 || m_types & 1 << i) && (this->*fun[i])())
				{
					m_stop = type(1 << i);
					return true;
				}
			}
			return false;
		}

		//stop
		bool stop_criteria::stop_step_maximum(void) const
		{
			return m_solver->m_step > m_solver->m_step_max;
		}
		bool stop_criteria::stop_load_limit_maximum(void) const
		{
			return m_solver->m_p_new > m_p_max;
		}
		bool stop_criteria::stop_load_limit_minimum(void) const
		{
			return m_solver->m_p_new < m_p_min;
		}
		bool stop_criteria::stop_load_local_minimum(void) const
		{
			const uint32_t s = m_solver->m_step;
			const double* p = m_solver->m_p_data;
			return s >= 3 && p[s - 2] < p[s - 3] && p[s - 2] < p[s - 1];
		}
		bool stop_criteria::stop_load_local_maximum(void) const
		{
			const uint32_t s = m_solver->m_step;
			const double* p = m_solver->m_p_data;
			return s >= 3 && p[s - 2] > p[s - 3] && p[s - 2] > p[s - 1];
		}
		bool stop_criteria::stop_load_value_negative(void) const
		{
			const uint32_t s = m_solver->m_step;
			const double* p = m_solver->m_p_data;
			return s >= 2 && p[s - 1] < 0 && p[s - 2] > 0;
		}
		bool stop_criteria::stop_load_value_positive(void) const
		{
			const uint32_t s = m_solver->m_step;
			const double* p = m_solver->m_p_data;
			return s >= 2 && p[s - 1] > 0 && p[s - 2] < 0;
		}
		bool stop_criteria::stop_state_limit_minimum(void) const
		{
			return m_solver->m_x_new[m_solver->m_watch_dof] < m_x_min;
		}
		bool stop_criteria::stop_state_limit_maximum(void) const
		{
			return m_solver->m_x_new[m_solver->m_watch_dof] > m_x_max;
		}
		bool stop_criteria::stop_state_local_minimum(void) const
		{
			const uint32_t s = m_solver->m_step;
			const uint32_t n = m_solver->m_size;
			const uint32_t w = m_solver->m_watch_dof;
			const double* x = m_solver->m_x_data + w;
			return s >= 3 && x[(s - 2) * n] < x[(s - 3) * n] && x[(s - 2) * n] < x[(s - 1) * n];
		}
		bool stop_criteria::stop_state_local_maximum(void) const
		{
			const uint32_t s = m_solver->m_step;
			const uint32_t n = m_solver->m_size;
			const uint32_t w = m_solver->m_watch_dof;
			const double* x = m_solver->m_x_data + w;
			return s >= 3 && x[(s - 2) * n] > x[(s - 3) * n] && x[(s - 2) * n] > x[(s - 1) * n];
		}
		bool stop_criteria::stop_state_value_negative(void) const
		{
			const uint32_t s = m_solver->m_step;
			const uint32_t n = m_solver->m_size;
			const uint32_t w = m_solver->m_watch_dof;
			const double* x = m_solver->m_x_data + w;
			return s >= 2 && x[(s - 1) * n] < 0 && x[(s - 2) * n] > 0;
		}
		bool stop_criteria::stop_state_value_positive(void) const
		{
			const uint32_t s = m_solver->m_step;
			const uint32_t n = m_solver->m_size;
			const uint32_t w = m_solver->m_watch_dof;
			const double* x = m_solver->m_x_data + w;
			return s >= 2 && x[(s - 1) * n] > 0 && x[(s - 2) * n] < 0;
		}
	}
}
