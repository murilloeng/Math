//Math
#include "Math/inc/Eigen/Dense.hpp"

namespace math
{
	namespace eigen
	{
		//constructor
		Dense::Dense(uint32_t order) : Base{order}
		{
			return;
		}
		
		//destructor
		Dense::~Dense(void)
		{
			return;
		}
	}
}