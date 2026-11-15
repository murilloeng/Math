#pragma once

//std
#include <cmath>
#include <cstdint>

//math
#include "MatrixFixed.hpp"

namespace math
{
	template<uint32_t n> class VectorFixed : public MatrixFixed<n, 1>
	{
	public:
		//constructor
		VectorFixed(void) = default;
		VectorFixed(const MatrixFixed<n, 1>&);

		explicit VectorFixed(double*);
		explicit VectorFixed(const double*);

		//destructor
		~VectorFixed(void) = default;

		//operators
		VectorFixed& operator=(const MatrixFixed<n, 1>&);

		//linear
		double norm(void) const;
		double inner(const VectorFixed&) const;
		VectorFixed cross(const VectorFixed&) const;
	};

	//constructor
	template<uint32_t n> VectorFixed<n>::VectorFixed(const MatrixFixed<n, 1>& M) : MatrixFixed<n, 1>(M)
	{
		return;
	}
	template<uint32_t n> VectorFixed<n>::VectorFixed(double* data_ptr) : MatrixFixed<n, 1>(data_ptr)
	{
		return;
	}
	template<uint32_t n> VectorFixed<n>::VectorFixed(const double* data_ref) : MatrixFixed<n, 1>(data_ref)
	{
		return;
	}

	//operators
	template<uint32_t n> VectorFixed<n>& VectorFixed<n>::operator=(const MatrixFixed<n, 1>& M)
	{
		MatrixFixed<n, 1>::operator=(M);
		return *this;
	}

	//linear
	template<uint32_t n> double VectorFixed<n>::norm(void) const
	{
		return sqrt(inner(*this));
	}
	template<uint32_t n> double VectorFixed<n>::inner(const VectorFixed& v) const
	{
		double s = 0;
		for(uint32_t i = 0; i < n; i++) s += this->m_data_ref[i] * v.m_data_ref[i];
		return s;
	}
	template<uint32_t n> VectorFixed<n> VectorFixed<n>::cross(const VectorFixed& v) const
	{
		VectorFixed w;
		static_assert(n == 3, "Error: Cross product is only defined for 3D vectors");
		w.m_data_ptr[0] = m_data_ref[1] * v.m_data_ref[2] - m_data_ref[2] * v.m_data_ref[1];
		w.m_data_ptr[1] = m_data_ref[2] * v.m_data_ref[0] - m_data_ref[0] * v.m_data_ref[2];
		w.m_data_ptr[2] = m_data_ref[0] * v.m_data_ref[1] - m_data_ref[1] * v.m_data_ref[0];
		return w;
	}
}