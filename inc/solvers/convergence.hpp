#pragma once

//std
#include <cstdint>

namespace math
{
	namespace solvers
	{
		class solver;
	}
}

namespace math
{
	namespace solvers
	{
		class convergence
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
			convergence(void);

			//destructor
			~convergence(void);

			//check
			bool check(void) const;

			//data
			type m_type;
			solver* m_solver;
			double m_tolerance;
		};
	}
}