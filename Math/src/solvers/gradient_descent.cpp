//std
#include <stdexcept>

//math
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/gradient_descent.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		gradient_descent::gradient_descent(void)
		{
			m_step_max = 1;
			m_attempt_max = 1;
			m_iteration_max = 100;
			m_convergence.m_type = convergence::type::fixed;
		}
		
		//destructor
		gradient_descent::~gradient_descent(void)
		{
			return;
		}

		//data
		uint32_t gradient_descent::state_set(void) const
		{
			return uint32_t(state::x);
		}
		uint32_t gradient_descent::force_set(void) const
		{
			return uint32_t(force::r);
		}
		uint32_t gradient_descent::tangent_set(void) const
		{
			return 0;
		}

		//solve
		void gradient_descent::check(void)
		{
			if(!m_gradient)
			{
				throw std::runtime_error("gradient descent called with gradient missing!");
			}
		}
		void gradient_descent::compute(void)
		{
			m_gradient(m_r, m_x_new);
		}
		void gradient_descent::predictor(void)
		{
			math::vector(m_dx, m_size).zeros();
		}
		void gradient_descent::corrector(void)
		{
			math::vector(m_dx, m_size) -= m_step_size * math::vector(m_r, m_size);
		}
	}
}