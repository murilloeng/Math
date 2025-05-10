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
			m_types(0), m_p_min(-DBL_MAX), m_p_max(+DBL_MAX), m_x_min(-DBL_MAX), m_x_max(+DBL_MAX)
		{
			return;
		}

		//destructor
		stop_criteria::~stop_criteria(void)
		{
			return;
		}

		//stop
		bool stop_criteria::stop(void) const
		{
			return 
				m_solver->m_p_new < m_p_min || 
				m_solver->m_p_new > m_p_max ||
				m_solver->m_step > m_solver->m_step_max ||
				m_solver->m_x_new[m_solver->m_watch_dof] < m_x_min ||
				m_solver->m_x_new[m_solver->m_watch_dof] > m_x_max ||
				(m_types & uint32_t(type::load_maximum) && stop_load_maximum()) ||
				(m_types & uint32_t(type::load_minimum) && stop_load_minimum()) ||
				(m_types & uint32_t(type::load_positive) && stop_load_positive()) ||
				(m_types & uint32_t(type::load_negative) && stop_load_negative()) ||
				(m_types & uint32_t(type::state_maximum) && stop_state_maximum()) ||
				(m_types & uint32_t(type::state_minimum) && stop_state_minimum()) ||
				(m_types & uint32_t(type::state_positive) && stop_state_positive()) ||
				(m_types & uint32_t(type::state_negative) && stop_state_negative());
		}

		//stop
		bool stop_criteria::stop_load_maximum(void) const
		{
			const uint32_t s = m_solver->m_step;
			const double* p = m_solver->m_p_data;
			return s >= 3 && p[s - 2] > p[s - 3] && p[s - 2] > p[s - 1];
		}
		bool stop_criteria::stop_load_minimum(void) const
		{
			const uint32_t s = m_solver->m_step;
			const double* p = m_solver->m_p_data;
			return s >= 3 && p[s - 2] < p[s - 3] && p[s - 2] < p[s - 1];
		}
		bool stop_criteria::stop_load_positive(void) const
		{
			const uint32_t s = m_solver->m_step;
			const double* p = m_solver->m_p_data;
			return s >= 2 && p[s - 1] > 0 && p[s - 2] < 0;
		}
		bool stop_criteria::stop_load_negative(void) const
		{
			const uint32_t s = m_solver->m_step;
			const double* p = m_solver->m_p_data;
			return s >= 2 && p[s - 1] < 0 && p[s - 2] > 0;
		}
		bool stop_criteria::stop_state_maximum(void) const
		{
			const uint32_t s = m_solver->m_step;
			const uint32_t n = m_solver->m_size;
			const uint32_t w = m_solver->m_watch_dof;
			const double* x = m_solver->m_x_data + w;
			return s >= 3 && x[(s - 2) * n] > x[(s - 3) * n] && x[(s - 2) * n] > x[(s - 1) * n];
		}
		bool stop_criteria::stop_state_minimum(void) const
		{
			const uint32_t s = m_solver->m_step;
			const uint32_t n = m_solver->m_size;
			const uint32_t w = m_solver->m_watch_dof;
			const double* x = m_solver->m_x_data + w;
			return s >= 3 && x[(s - 2) * n] < x[(s - 3) * n] && x[(s - 2) * n] < x[(s - 1) * n];
		}
		bool stop_criteria::stop_state_positive(void) const
		{
			const uint32_t s = m_solver->m_step;
			const uint32_t n = m_solver->m_size;
			const uint32_t w = m_solver->m_watch_dof;
			const double* x = m_solver->m_x_data + w;
			return s >= 2 && x[(s - 1) * n] > 0 && x[(s - 2) * n] < 0;
		}
		bool stop_criteria::stop_state_negative(void) const
		{
			const uint32_t s = m_solver->m_step;
			const uint32_t n = m_solver->m_size;
			const uint32_t w = m_solver->m_watch_dof;
			const double* x = m_solver->m_x_data + w;
			return s >= 2 && x[(s - 1) * n] < 0 && x[(s - 2) * n] > 0;
		}
	}
}
