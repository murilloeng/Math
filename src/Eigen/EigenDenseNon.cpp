//std
#include <cfloat>

//Math
#include "Math/inc/Eigen/EigenDenseNon.hpp"

namespace math
{
	namespace eigen
	{
		//constructor
		EigenDenseNon::EigenDenseNon(uint32_t order) : EigenDense{order}
		{
			return;
		}
		
		//destructor
		EigenDenseNon::~EigenDenseNon(void)
		{
			return;
		}
	}
}