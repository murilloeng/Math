#pragma once

//std
#include <cstdint>

namespace math
{
	namespace solvers
	{
		class Solver;
	}
}

namespace math
{
	namespace solvers
	{
		class Convergence
		{
		public:
			//types
			enum class type : uint32_t
			{
				fixed = 1 << 0, 
				force = 1 << 1,
				last
			};

			//constructor
			Convergence(void);

			//destructor
			~Convergence(void);

			//check
			bool check(void) const;

			//data
			type m_type;
			Solver* m_solver;
			double m_tolerance;
		};
	}
}