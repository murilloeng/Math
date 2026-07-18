//Math
#include "Math/inc/Eigen/EigenDense.hpp"

namespace math
{
	namespace eigen
	{
		//constructor
		EigenDense::EigenDense(uint32_t order) : Eigen{order}
		{
			return;
		}
		
		//destructor
		EigenDense::~EigenDense(void)
		{
			return;
		}
	}
}