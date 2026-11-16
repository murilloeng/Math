#pragma once

//std
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <initializer_list>

//Math
#include "Math/inc/Linear/Matrix.hpp"

//defines
#ifndef VECTOR_STATIC_SIZE
#define VECTOR_STATIC_SIZE 18
#endif

namespace math
{
	class SpanN;
	class MatrixN;
}

namespace math
{
	class VectorN
	{
	public:
		//constructors
		VectorN(const VectorN&);
		explicit VectorN(uint32_t);
		VectorN(double*, uint32_t);
		VectorN(const double*, uint32_t);
		VectorN(std::initializer_list<double>);
		template<uint32_t n> VectorN(const Vector<n>&);

		//destructor
		virtual ~VectorN(void);

		//size
		uint32_t rows(void) const;

		//data
		double* data(void);
		const double* data(void) const;

		//direction
		double norm(void) const;
		VectorN& normalize(void);
		VectorN unit(void) const;

		//utility
		VectorN& zeros(void);
		VectorN& randu(double = - 1, double = +1);

		//print
		void print(const char* = "", double = 0) const;

		//products
		MatrixN outer(void) const;
		MatrixN outer(const VectorN&) const;
		double inner(const VectorN&) const;
		double inner(const double*) const;

		//operators
		VectorN operator+(void) const;
		VectorN operator-(void) const;
		VectorN operator/(double) const;

		VectorN operator+(const VectorN&) const;
		VectorN operator-(const VectorN&) const;

		VectorN& operator=(double);
		VectorN& operator=(const double*);
		VectorN& operator=(const VectorN&);
		template<uint32_t n> VectorN& operator=(const Vector<n>&);
		VectorN& operator=(std::initializer_list<double>);

		VectorN& operator+=(double);
		VectorN& operator-=(double);
		VectorN& operator*=(double);
		VectorN& operator/=(double);

		VectorN& operator+=(const double*);
		VectorN& operator-=(const double*);
		VectorN& operator+=(const VectorN&);
		VectorN& operator-=(const VectorN&);
		template<uint32_t n> VectorN& operator+=(const Vector<n>&);
		template<uint32_t n> VectorN& operator-=(const Vector<n>&);

		double& operator[](uint32_t);
		const double& operator[](uint32_t) const;

		friend VectorN operator*(double, const VectorN&);

	private:
		//data
		bool m_own;
		uint32_t m_rows;
		double* m_data_ptr;
		const double* m_data_ref;
		double m_data_mem[VECTOR_STATIC_SIZE];

		//friends
		friend class math::SpanN;
		friend class math::MatrixN;
	};
}

namespace math
{
	//constructors
	template<uint32_t n> VectorN::VectorN(const Vector<n>& v) : VectorN(n)
	{
		memcpy(m_data_ptr, v.data(), n * sizeof(double));
	}

	//operators
	template<uint32_t n> VectorN& VectorN::operator=(const Vector<n>& v)
	{
		//check
		if(m_rows != n)
		{
			throw std::runtime_error("VectorN assign has incompatible dimensions!");
		}
		//assign
		memcpy(m_data_ptr, v.data(), n * sizeof(double));
		//return
		return *this;
	}
	template<uint32_t n> VectorN& VectorN::operator+=(const Vector<n>& v)
	{
		if(m_rows != n)
		{
			throw std::runtime_error("VectorN increment has incompatible dimensions!");
		}
		for(uint32_t i = 0; i < n; i++)
		{
			m_data_ptr[i] += v[i];
		}
		return *this;
	}
	template<uint32_t n> VectorN& VectorN::operator-=(const Vector<n>& v)
	{
		if(m_rows != n)
		{
			throw std::runtime_error("VectorN increment has incompatible dimensions!");
		}
		for(uint32_t i = 0; i < n; i++)
		{
			m_data_ptr[i] -= v[i];
		}
		return *this;
	}
}
