//Math
#include "Math/inc/Eigen/EigenDenseStd.hpp"

namespace math
{
	namespace eigen
	{
		//constructor
		EigenDenseStd::EigenDenseStd(uint32_t order, double* A) : EigenDense{order}, m_A{A}
		{
			return;
		}
		
		//destructor
		EigenDenseStd::~EigenDenseStd(void)
		{
			return;
		}
	}
}