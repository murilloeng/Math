//Math
#include "Math/inc/Eigen/Base.hpp"

namespace math
{
	namespace eigen
	{
		//constructor
		Base::Base(uint32_t order) : m_order{order}
		{
			return;
		}
		
		//destructor
		Base::~Base(void)
		{
			return;
		}
	}
}