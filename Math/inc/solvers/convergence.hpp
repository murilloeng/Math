#pragma once

//std
#include <cstdint>

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
			double* m_r;
			double* m_g;
			uint32_t m_size;
			double m_tolerance;
		};
	}
}