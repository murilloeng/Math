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