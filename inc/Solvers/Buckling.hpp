#pragma once

//std


//Math
#include "Math/inc/Solvers/Solver.hpp"

namespace math
{
	namespace solvers
	{
		class Buckling : public Solver
		{
		public:
			//constructor
			Buckling(void);

			//destructor
			~Buckling(void);

			//data
			uint32_t state_set(void) const override;
			uint32_t force_set(void) const override;
			uint32_t tangent_set(void) const override;

			//solve
			void solve(void) override;

			//data
			bool m_full;
			uint32_t m_modes;
			std::function<void(double*, const double*)> m_stiffness;
		};
	}
}