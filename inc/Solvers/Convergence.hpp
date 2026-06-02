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
			enum class Type : uint32_t
			{
				Fixed = 1 << 0, 
				Force = 1 << 1,
				Last
			};

			//constructor
			Convergence(void);

			//destructor
			~Convergence(void);

			//check
			bool check(void) const;

			//data
			Type m_type;
			Solver* m_solver;
			double m_tolerance;
		};
	}
}