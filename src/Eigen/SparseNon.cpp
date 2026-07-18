//Math
#include "Math/inc/Eigen/SparseNon.hpp"

namespace math
{
	namespace eigen
	{
		//constructor
		SparseNon::SparseNon(uint32_t order, uint32_t modes, uint32_t vectors, const int32_t* rows_map, const int32_t* cols_map) : 
			Sparse{order, modes, vectors, rows_map, cols_map}
		{
			return;
		}
		
		//destructor
		SparseNon::~SparseNon(void)
		{
			return;
		}
	}
}