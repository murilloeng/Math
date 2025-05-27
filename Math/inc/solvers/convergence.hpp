#pragma once

//std
#include <cstdint>

namespace math
{
	namespace solvers
	{
		class newton_raphson;
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
			convergence(const newton_raphson*);

			//destructor
			~convergence(void);

			//check
			bool check(void) const;

			//data
			type m_type;
			double m_tolerance;
			const newton_raphson* m_solver;
		};
	}
}