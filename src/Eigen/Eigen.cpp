//Math
#include "Math/inc/Eigen/Eigen.hpp"

namespace math
{
	namespace eigen
	{
		//constructor
		Eigen::Eigen(uint32_t order) : m_order{order}
		{
			return;
		}
		
		//destructor
		Eigen::~Eigen(void)
		{
			return;
		}
	}
}