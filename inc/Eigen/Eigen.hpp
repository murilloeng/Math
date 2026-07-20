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

			//data
			uint32_t order(uint32_t);
			uint32_t order(void) const;

			double tolerance(double);
			double tolerance(void) const;

			//compute
			virtual bool compute(void) = 0;

		protected:
			//data
			uint32_t m_order;
			double m_tolerance;
		};
	}
}