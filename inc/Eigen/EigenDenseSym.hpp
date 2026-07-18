#pragma once

//Math
#include "Math/inc/Eigen/EigenDense.hpp"

namespace math
{
	namespace eigen
	{
		class EigenDenseSym : public EigenDense
		{
		public:
			//constructor
			EigenDenseSym(uint32_t);

			//destructor
			virtual ~EigenDenseSym(void);

			//data
			char range(char);
			char range(void) const;

			uint32_t modes(void) const;

			double value_min(double);
			double value_min(void) const;

			double value_max(double);
			double value_max(void) const;

			uint32_t index_min(uint32_t);
			uint32_t index_min(void) const;

			uint32_t index_max(uint32_t);
			uint32_t index_max(void) const;

		protected:
			//data
			char m_range;
			uint32_t m_modes;
			double m_value_min;
			double m_value_max;
			uint32_t m_index_min;
			uint32_t m_index_max;
		};
	}
}