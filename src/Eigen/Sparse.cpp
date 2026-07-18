//Math
#include "Math/inc/Eigen/Sparse.hpp"

namespace math
{
	namespace eigen
	{
		//constructor
		Sparse::Sparse(uint32_t order, uint32_t modes, uint32_t vectors, const int32_t* rows_map, const int32_t* cols_map) : 
			Base{order}, m_modes{modes}, m_vectors{vectors}, m_rows_map{rows_map}, m_cols_map{cols_map}
		{
			return;
		}
		
		//destructor
		Sparse::~Sparse(void)
		{
			return;
		}

		//data
		uint32_t Sparse::modes(void) const
		{
			return m_modes;
		}
		uint32_t Sparse::modes(uint32_t modes)
		{
			return m_modes = modes;
		}

		uint32_t Sparse::vectors(void) const
		{
			return m_vectors;
		}
		uint32_t Sparse::vectors(uint32_t vectors)
		{
			return m_vectors = vectors;
		}

		const int32_t* Sparse::rows_map(void) const
		{
			return m_rows_map;
		}
		const int32_t* Sparse::rows_map(const int32_t* rows_map)
		{
			return m_rows_map = rows_map;
		}

		const int32_t* Sparse::cols_map(void) const
		{
			return m_cols_map;
		}
		const int32_t* Sparse::cols_map(const int32_t* cols_map)
		{
			return m_cols_map = cols_map;
		}
	}
}