#pragma once

//Math
#include "Math/inc/Solvers/Implicit.hpp"

namespace math
{
	namespace solvers
	{
		class GradientDescent : virtual public Implicit
		{
		public:
			//constructor
			GradientDescent(void);

			//destructor
			~GradientDescent(void);

			//types
			typedef std::function<void(double*, const double*)> Gradient;

			//solve
			void solve(void) override;

			//data
			double step_size(double);
			double step_size(void) const;

			Gradient gradient(Gradient);
			Gradient gradient(void) const;
			void gradient(double*, const double*) const;

			//data
			uint32_t state_set(void) const override;
			uint32_t force_set(void) const override;
			uint32_t tangent_set(void) const override;

		private:
			//solve
			void check(void) override;
			void compute(void) override;
			void predictor(void) override;
			void corrector(void) override;

			//data
			using Implicit::m_attempt, Implicit::m_attempt_max;
			using Implicit::m_convergence, Implicit::m_continuation;

			//data
			double m_step_size;
			Gradient m_gradient;
		};
	}
}