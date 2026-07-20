//std
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Eigen/Sparse.hpp"

namespace math
{
	namespace eigen
	{
		//constructor
		Sparse::Sparse(uint32_t order, uint32_t modes, uint32_t vectors, const int32_t* rows_map, const int32_t* cols_map) : 
			Eigen{order}, m_mode{1U}, m_modes{modes}, m_vectors{vectors}, m_iteration_max{300U}, 
			m_type{"SA"}, m_shift{0}, m_rows_map{rows_map}, m_cols_map{cols_map}
		{
			return;
		}
		
		//destructor
		Sparse::~Sparse(void)
		{
			return;
		}

		//data
		uint32_t Sparse::mode(void) const
		{
			return m_mode;
		}
		uint32_t Sparse::mode(uint32_t type)
		{
			return m_mode = type;
		}

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

		uint32_t Sparse::iteration_max(void) const
		{
			return m_iteration_max;
		}
		uint32_t Sparse::iteration_max(uint32_t iteration_max)
		{
			return m_iteration_max = iteration_max;
		}

		double Sparse::shift(void) const
		{
			return m_shift;
		}
		double Sparse::shift(double shift)
		{
			return m_shift = shift;
		}

		const char* Sparse::type(void) const
		{
			return m_type;
		}
		const char* Sparse::type(const char* type)
		{
			if(strcmp(type, "LA") && strcmp(type, "SA") && strcmp(type, "LM") && strcmp(type, "SM") && strcmp(type, "BE"))
			{
				throw std::runtime_error("Error: Eigen Sparse type not avaliable!");
			}
			return strcpy(m_type, type);
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