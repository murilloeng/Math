//std
#include <cfloat>

//Math
#include "Math/inc/Eigen/EigenDenseSym.hpp"

namespace math
{
	namespace eigen
	{
		//constructor
		EigenDenseSym::EigenDenseSym(uint32_t order) : EigenDense{order},
			m_range{'A'}, m_value_min{-DBL_MAX}, m_value_max{+DBL_MAX}, m_index_min{0}, m_index_max{order}
		{
			return;
		}
		
		//destructor
		EigenDenseSym::~EigenDenseSym(void)
		{
			return;
		}

		//data
		char EigenDenseSym::range(void) const
		{
			return m_range;
		}
		char EigenDenseSym::range(char range)
		{
			return m_range = range;
		}

		uint32_t EigenDenseSym::modes(void) const
		{
			return m_modes;
		}

		double EigenDenseSym::value_min(void) const
		{
			return m_value_min;
		}
		double EigenDenseSym::value_min(double value_min)
		{
			return m_value_min = value_min;
		}

		double EigenDenseSym::value_max(void) const
		{
			return m_value_max;
		}
		double EigenDenseSym::value_max(double value_max)
		{
			return m_value_max = value_max;
		}

		uint32_t EigenDenseSym::index_min(void) const
		{
			return m_index_min;
		}
		uint32_t EigenDenseSym::index_min(uint32_t index_min)
		{
			return m_index_min = index_min;
		}

		uint32_t EigenDenseSym::index_max(void) const
		{
			return m_index_max;
		}
		uint32_t EigenDenseSym::index_max(uint32_t index_max)
		{
			return m_index_max = index_max;
		}
	}
}