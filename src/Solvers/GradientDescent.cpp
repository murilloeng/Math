//std
#include <stdexcept>

//Math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Solvers/GradientDescent.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		GradientDescent::GradientDescent(void)
		{
			// m_convergence.m_type = Convergence::Type::Fixed;
		}

		//destructor
		GradientDescent::~GradientDescent(void)
		{
			return;
		}

		//data
		uint32_t GradientDescent::state_set(void) const
		{
			return uint32_t(State::x);
		}
		uint32_t GradientDescent::force_set(void) const
		{
			return uint32_t(Force::r);
		}
		uint32_t GradientDescent::tangent_set(void) const
		{
			return 0;
		}

		//solve
		void GradientDescent::check(void)
		{
			if(!m_gradient)
			{
				throw std::runtime_error("gradient descent called with gradient missing!");
			}
		}
		void GradientDescent::compute(void)
		{
			m_gradient(m_r, m_x_new);
		}
		void GradientDescent::predictor(void)
		{
			math::Vector(m_dx, m_size).zeros();
		}
		void GradientDescent::corrector(void)
		{
			math::Vector(m_dx, m_size) -= m_step_size * math::Vector(m_r, m_size);
		}
	}
}