//std
#include <cfloat>

//Math
#include "Math/inc/Solvers/Incremental.hpp"
#include "Math/inc/Solvers/StopCriteria.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		StopCriteria::StopCriteria(Incremental* solver) :
			m_stop{Type::Last}, m_types{uint32_t(Type::StepMaximum) & uint32_t(Type::TimeMaximum)}, m_solver{solver},
			m_state_min{-DBL_MAX}, m_parameter_min{-DBL_MAX}, m_state_max{+DBL_MAX}, m_parameter_max{+DBL_MAX}
		{
			return;
		}

		//destructor
		StopCriteria::~StopCriteria(void)
		{
			return;
		}

		//stop
		bool StopCriteria::stop(void)
		{
			//data
			bool(StopCriteria::*fun[])(void) const = {
				&StopCriteria::stop_step_maximum,
				&StopCriteria::stop_time_maximum,
				&StopCriteria::stop_state_limit_minimum, &StopCriteria::stop_state_limit_maximum,
				&StopCriteria::stop_state_local_minimum, &StopCriteria::stop_state_local_maximum,
				&StopCriteria::stop_state_value_negative, &StopCriteria::stop_state_value_positive,
				&StopCriteria::stop_parameter_limit_minimum, &StopCriteria::stop_parameter_limit_maximum,
				&StopCriteria::stop_parameter_local_minimum, &StopCriteria::stop_parameter_local_maximum,
				&StopCriteria::stop_parameter_value_negative, &StopCriteria::stop_parameter_value_positive
			};
			//stop
			m_stop = Type::Last;
			for(uint32_t i = 0; i < uint32_t(Type::Last); i++)
			{
				if((i < 2 || m_types & 1 << i) && (this->*fun[i])())
				{
					m_stop = Type(1 << i);
					return true;
				}
			}
			return false;
		}

		//data
		uint32_t StopCriteria::types(void) const
		{
			return m_types;
		}
		uint32_t StopCriteria::add_type(Type type)
		{
			return m_types |= 1 << uint32_t(type);
		}

		double StopCriteria::state_min(void) const
		{
			return m_state_min;
		}
		double StopCriteria::state_min(double state_min)
		{
			return m_state_min = state_min;
		}

		double StopCriteria::state_max(void) const
		{
			return m_state_max;
		}
		double StopCriteria::state_max(double state_max)
		{
			return m_state_max = state_max;
		}

		double StopCriteria::parameter_min(void) const
		{
			return m_parameter_min;
		}
		double StopCriteria::parameter_min(double parameter_min)
		{
			return m_parameter_min = parameter_min;
		}

		double StopCriteria::parameter_max(void) const
		{
			return m_parameter_max;
		}
		double StopCriteria::parameter_max(double parameter_max)
		{
			return m_parameter_max = parameter_max;
		}

		//stop
		bool StopCriteria::stop_step_maximum(void) const
		{
			return m_solver->m_step > m_solver->m_step_max;
		}
		bool StopCriteria::stop_time_maximum(void) const
		{
			return m_solver->m_t_new > m_solver->m_t_max;
		}
		bool StopCriteria::stop_state_limit_minimum(void) const
		{
			return m_solver->m_x_new[m_solver->m_watch_dof] < m_state_min;
		}
		bool StopCriteria::stop_state_limit_maximum(void) const
		{
			return m_solver->m_x_new[m_solver->m_watch_dof] > m_state_max;
		}
		bool StopCriteria::stop_state_local_minimum(void) const
		{
			const uint32_t s = m_solver->m_step;
			const uint32_t n = m_solver->m_size;
			const uint32_t w = m_solver->m_watch_dof;
			const double* x = m_solver->m_x_data + w;
			return s >= 3 && x[(s - 2) * n] < x[(s - 3) * n] && x[(s - 2) * n] < x[(s - 1) * n];
		}
		bool StopCriteria::stop_state_local_maximum(void) const
		{
			const uint32_t s = m_solver->m_step;
			const uint32_t n = m_solver->m_size;
			const uint32_t w = m_solver->m_watch_dof;
			const double* x = m_solver->m_x_data + w;
			return s >= 3 && x[(s - 2) * n] > x[(s - 3) * n] && x[(s - 2) * n] > x[(s - 1) * n];
		}
		bool StopCriteria::stop_state_value_negative(void) const
		{
			const uint32_t s = m_solver->m_step;
			const uint32_t n = m_solver->m_size;
			const uint32_t w = m_solver->m_watch_dof;
			const double* x = m_solver->m_x_data + w;
			return s >= 2 && x[(s - 1) * n] < 0 && x[(s - 2) * n] > 0;
		}
		bool StopCriteria::stop_state_value_positive(void) const
		{
			const uint32_t s = m_solver->m_step;
			const uint32_t n = m_solver->m_size;
			const uint32_t w = m_solver->m_watch_dof;
			const double* x = m_solver->m_x_data + w;
			return s >= 2 && x[(s - 1) * n] > 0 && x[(s - 2) * n] < 0;
		}
		bool StopCriteria::stop_parameter_limit_minimum(void) const
		{
			return m_solver->m_p_new < m_parameter_min;
		}
		bool StopCriteria::stop_parameter_limit_maximum(void) const
		{
			return m_solver->m_p_new > m_parameter_max;
		}
		bool StopCriteria::stop_parameter_local_minimum(void) const
		{
			const uint32_t s = m_solver->m_step;
			const double* p = m_solver->m_p_data;
			return s >= 3 && p[s - 2] < p[s - 3] && p[s - 2] < p[s - 1];
		}
		bool StopCriteria::stop_parameter_local_maximum(void) const
		{
			const uint32_t s = m_solver->m_step;
			const double* p = m_solver->m_p_data;
			return s >= 3 && p[s - 2] > p[s - 3] && p[s - 2] > p[s - 1];
		}
		bool StopCriteria::stop_parameter_value_negative(void) const
		{
			const uint32_t s = m_solver->m_step;
			const double* p = m_solver->m_p_data;
			return s >= 2 && p[s - 1] < 0 && p[s - 2] > 0;
		}
		bool StopCriteria::stop_parameter_value_positive(void) const
		{
			const uint32_t s = m_solver->m_step;
			const double* p = m_solver->m_p_data;
			return s >= 2 && p[s - 1] > 0 && p[s - 2] < 0;
		}
	}
}
