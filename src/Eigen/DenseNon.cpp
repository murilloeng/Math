//std
#include <cfloat>

//Math
#include "Math/inc/Eigen/DenseNon.hpp"

namespace math
{
	namespace eigen
	{
		//constructor
		DenseNon::DenseNon(uint32_t order) : Dense{order}
		{
			return;
		}
		
		//destructor
		DenseNon::~DenseNon(void)
		{
			return;
		}
	}
}