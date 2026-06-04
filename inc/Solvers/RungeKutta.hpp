#pragma once

//Math
#include "Math/inc/Solvers/Solver.hpp"

namespace math
{
	namespace solvers
	{
		class RungeKutta : public Solver
		{
		public:
			//constructors
			RungeKutta(void);

			//destructor
			~RungeKutta(void);

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

		public:
			//data
			void (*m_system_1)(double*, const double*, double);
			void (*m_system_2)(double*, const double*, const double*, double);
		};
	}
}