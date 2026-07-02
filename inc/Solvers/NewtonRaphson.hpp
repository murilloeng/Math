#pragma once

//std
#include <cstdint>
#include <functional>

//Math
#include "Math/inc/Solvers/Implicit.hpp"
#include "Math/inc/Solvers/Incremental.hpp"

namespace math
{
	namespace solvers
	{
		class NewtonRaphson : virtual public Implicit, virtual public Incremental
		{
		public:
			//constructors
			NewtonRaphson(void);

			//destructor
			~NewtonRaphson(void);

			//data
			uint32_t state_set(void) const override;
			uint32_t force_set(void) const override;
			uint32_t tangent_set(void) const override;

		protected:
			//solve
			void check(void) override;
			void print(void) override;
			void setup(void) override;
			void compute(void) override;
			void predictor(void) override;
			void corrector(void) override;

		public:
			//data
			std::function<void(double*, double, const double*)> m_residue;
			std::function<void(double*, double, const double*)> m_tangent_1;
			std::function<void(double*, double, const double*)> m_tangent_2;
			std::function<void(double*, double*, const double*)> m_system_1;
			std::function<void(double*, double*, double*, double, const double*)> m_system_2;
		};
	}
}