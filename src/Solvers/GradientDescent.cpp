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
			m_convergence.m_type = Convergence::Type::Fixed;
		}

		//destructor
		GradientDescent::~GradientDescent(void)
		{
			return;
		}

		//solve
		void GradientDescent::solve(void)
		{
			check();
			setup();
			predictor();
			for(m_iteration = 0; m_iteration < m_iteration_max; m_iteration++)
			{
				apply();
				compute();
				if(equilibrium()) break; else corrector();
			}
		}

		//data
		double GradientDescent::step_size(void) const
		{
			return m_step_size;
		}
		double GradientDescent::step_size(double step_size)
		{
			return m_step_size = step_size;
		}

		GradientDescent::Gradient GradientDescent::gradient(void) const
		{
			return m_gradient;
		}
		GradientDescent::Gradient GradientDescent::gradient(Gradient gradient)
		{
			return m_gradient = gradient;
		}
		void GradientDescent::gradient(double* g, const double* x) const
		{
			m_gradient(g, x);
		}

		//data
		uint32_t GradientDescent::state_set(void) const
		{
			return 1 << uint32_t(State::x);
		}
		uint32_t GradientDescent::force_set(void) const
		{
			return 1 << uint32_t(Force::r);
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