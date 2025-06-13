#pragma once

//std
#include <cstdint>
#include <functional>

//math
#include "Math/Math/inc/solvers/solver.hpp"

//x: state vector
//r: residual force vector
//p: continuation parameter

//target system: r(x, p) = 0

//tangent on x: K(x, p) = -dr/dx(x, p)
//tangent on p: fe(x, p) = +dr/dp(x, p)

namespace math
{
	namespace solvers
	{
		class newton_raphson : public solver
		{
		public:
			//constructors
			newton_raphson(void);

			//destructor
			~newton_raphson(void);

			//data
			void save(const char*) const;
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
			std::function<void(double*, double, const double*)> m_residue;
			std::function<void(double*, double, const double*)> m_tangent_1;
			std::function<void(double*, double, const double*)> m_tangent_2;
			std::function<void(double*, double*, const double*)> m_system_1;
			std::function<void(double*, double*, double*, double, const double*)> m_system_2;
		};
	}
}