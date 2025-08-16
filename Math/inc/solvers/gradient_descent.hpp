#pragma once

//math
#include "Math/Math/inc/solvers/solver.hpp"

namespace math
{
	namespace solvers
	{
		class gradient_descent : public solver
		{
		public:
			//constructor
			gradient_descent(void);

			//destructor
			~gradient_descent(void);

			//data
			uint32_t state_set(void) const override;
			uint32_t force_set(void) const override;
			uint32_t tangent_set(void) const override;

			//solve
			void check(void) override;
			void compute(void) override;
			void predictor(void) override;
			void corrector(void) override;

			//data
			double m_step_size;
			std::function<void(double*, const double*)> m_gradient;
		};
	}
}