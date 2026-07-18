#pragma once

//std
#include <cstdint>

namespace math
{
	namespace eigen
	{
		class Base
		{
		public:
			//constructor
			Base(uint32_t);

			//destructor
			virtual ~Base(void);

			//compute
			virtual bool compute(void) = 0;

		protected:
			//data
			uint32_t m_order;
		};
	}
}