#pragma once

//std
#include <cstdint>

namespace math
{
	namespace eigen
	{
		class Eigen
		{
		public:
			//constructor
			Eigen(uint32_t);

			//destructor
			virtual ~Eigen(void);

			//compute
			virtual bool compute(void) = 0;

		protected:
			//data
			uint32_t m_order;
		};
	}
}