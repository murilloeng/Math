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
			bool equilibrium(void) override;

			//compute
			void compute_tangent_1(void);
			void compute_tangent_2(void);
			void compute_tangent_3(void);
			void compute_tangent_4(void);

		public:
			//data
			std::function<void(double*, const double*)> m_inertia;
			std::function<void(double*, const double*, const double*)> m_internal_force;
			std::function<void(double*, const double*, const double*, double)> m_external_force;

		};
	}
}