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

			//solve
			void compute_step(void) override;

			//data
			uint32_t state_set(void) const override;
			uint32_t force_set(void) const override;
			uint32_t tangent_set(void) const override;

			//types
			typedef std::function<void(double*, double, const double*)> Function;
			typedef std::function<void(double*, double*, double*, double, const double*)> System;

			//data
			System system(System);
			System system(void) const;

			Function residue(Function);
			Function residue(void) const;

			Function tangent_p(Function);
			Function tangent_p(void) const;

			Function tangent_x(Function);
			Function tangent_x(void) const;

		protected:
			//solve
			void check(void) override;
			void print(void) override;
			void setup(void) override;
			void compute(void) override;
			void predictor(void) override;
			void corrector(void) override;

			//data
			System m_system;
			Function m_residue;
			Function m_tangent_p;
			Function m_tangent_x;
		};
	}
}