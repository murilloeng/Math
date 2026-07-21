//std
#include <cstring>

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
			strcpy(m_type, "SR");
		}
		
		//destructor
		SparseNon::~SparseNon(void)
		{
			return;
		}

		//data
		double SparseNon::shift(uint32_t index) const
		{
			return m_shift[index];
		}
		double SparseNon::shift(uint32_t index, double shift)
		{
			return m_shift[index] = shift;
		}
	}
}