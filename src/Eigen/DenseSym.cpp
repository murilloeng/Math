//std
#include <cfloat>

//Math
#include "Math/inc/Eigen/DenseSym.hpp"

namespace math
{
	namespace eigen
	{
		//constructor
		DenseSym::DenseSym(uint32_t order) : Dense{order},
			m_range{'A'}, m_value_min{-DBL_MAX}, m_value_max{+DBL_MAX}, m_index_min{0}, m_index_max{order}
		{
			return;
		}
		
		//destructor
		DenseSym::~DenseSym(void)
		{
			return;
		}

		//data
		char DenseSym::range(void) const
		{
			return m_range;
		}
		char DenseSym::range(char range)
		{
			return m_range = range;
		}

		uint32_t DenseSym::modes(void) const
		{
			return m_modes;
		}

		double DenseSym::value_min(void) const
		{
			return m_value_min;
		}
		double DenseSym::value_min(double value_min)
		{
			return m_value_min = value_min;
		}

		double DenseSym::value_max(void) const
		{
			return m_value_max;
		}
		double DenseSym::value_max(double value_max)
		{
			return m_value_max = value_max;
		}

		uint32_t DenseSym::index_min(void) const
		{
			return m_index_min;
		}
		uint32_t DenseSym::index_min(uint32_t index_min)
		{
			return m_index_min = index_min;
		}

		uint32_t DenseSym::index_max(void) const
		{
			return m_index_max;
		}
		uint32_t DenseSym::index_max(uint32_t index_max)
		{
			return m_index_max = index_max;
		}
	}
}