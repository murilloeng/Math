#pragma once

//std
#include <cstdint>

//Math
#include "Math/inc/Solvers/Implicit.hpp"
#include "Math/inc/Solvers/Incremental.hpp"

namespace math
{
	namespace solvers
	{
		class Newmark : virtual public Implicit, virtual public Incremental
		{
		public:
			//constructors
			Newmark(void);

			//destructor
			~Newmark(void);

			//data
			uint32_t state_set(void) const override;
			uint32_t force_set(void) const override;
			uint32_t tangent_set(void) const override;

		private:
			//solve
			void check(void) override;
			void print(void) override;
			void setup(void) override;
			void compute(void) override;
			void predictor(void) override;
			void corrector(void) override;

		public:
			//data
			double m_g, m_b;

			std::function<void(double*, const double*, const double*)> m_internal_force;
			std::function<void(double*, const double*, const double*, double)> m_external_force;

			std::function<void(double*, const double*)> m_inertia;
			std::function<void(double*, const double*, const double*, double)> m_damping;
			std::function<void(double*, const double*, const double*, const double*, double)> m_stiffness;
		};
	}
}