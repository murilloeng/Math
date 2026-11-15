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
		using MatrixFixed<n, 1>::MatrixFixed;
		VectorFixed(const MatrixFixed<n, 1>&);

		//destructor
		~VectorFixed(void) = default;

		//operators
		VectorFixed& operator=(const MatrixFixed<n, 1>&);

		//linear
		double norm(void) const;
		double inner(const VectorFixed&) const;
		double cross(const VectorFixed&) const requires(n == 2);
		VectorFixed cross(const VectorFixed&) const requires(n == 3);
	};

	//constructors
	template<uint32_t n> VectorFixed<n>::VectorFixed(const MatrixFixed<n, 1>& M) : MatrixFixed<n, 1>(M)
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
		for(uint32_t i = 0; i < n; i++) s += this->m_data[i] * v.m_data[i];
		return s;
	}
	template<uint32_t n> double VectorFixed<n>::cross(const VectorFixed& v) const requires(n == 2)
	{
		return this->m_data[0] * v.m_data[1] - this->m_data[1] * v.m_data[0];
	}
	template<uint32_t n> VectorFixed<n> VectorFixed<n>::cross(const VectorFixed& v) const requires(n == 3)
	{
		VectorFixed w;
		w.m_data[0] = this->m_data[1] * v.m_data[2] - this->m_data[2] * v.m_data[1];
		w.m_data[1] = this->m_data[2] * v.m_data[0] - this->m_data[0] * v.m_data[2];
		w.m_data[2] = this->m_data[0] * v.m_data[1] - this->m_data[1] * v.m_data[0];
		return w;
	}
}