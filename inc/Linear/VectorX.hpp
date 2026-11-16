#pragma once

//std
#include <cstdint>
#include <initializer_list>

//defines
#ifndef VECTOR_STATIC_SIZE
#define VECTOR_STATIC_SIZE 18
#endif

namespace math
{
	class SpanX;
	class MatrixX;
}

namespace math
{
	class VectorX
	{
	public:
		//constructors
		VectorX(const VectorX&);
		explicit VectorX(uint32_t);
		VectorX(double*, uint32_t);
		VectorX(const double*, uint32_t);
		VectorX(std::initializer_list<double>);

		//destructor
		virtual ~VectorX(void);

		//size
		uint32_t rows(void) const;

		//data
		double* data(void);
		const double* data(void) const;

		//direction
		double norm(void) const;
		VectorX& normalize(void);
		VectorX unit(void) const;
		
		//utility
		VectorX& zeros(void);
		VectorX& randu(double = - 1, double = +1);

		//print
		void print(const char* = "", double = 0) const;

		//products
		MatrixX outer(void) const;
		MatrixX outer(const VectorX&) const;
		double inner(const VectorX&) const;
		double inner(const double*) const;

		//operators
		VectorX operator+(void) const;
		VectorX operator-(void) const;
		VectorX operator/(double) const;

		VectorX operator+(const VectorX&) const;
		VectorX operator-(const VectorX&) const;

		VectorX& operator=(double);
		VectorX& operator=(const double*);
		VectorX& operator=(const VectorX&);
		VectorX& operator=(std::initializer_list<double>);

		VectorX& operator+=(double);
		VectorX& operator-=(double);
		VectorX& operator*=(double);
		VectorX& operator/=(double);

		VectorX& operator+=(const double*);
		VectorX& operator-=(const double*);
		VectorX& operator+=(const VectorX&);
		VectorX& operator-=(const VectorX&);

		double& operator[](uint32_t);
		const double& operator[](uint32_t) const;

		friend VectorX operator*(double, const VectorX&);

	private:
		//data
		bool m_own;
		uint32_t m_rows;
		double* m_data_ptr;
		const double* m_data_ref;
		double m_data_mem[VECTOR_STATIC_SIZE];

		//friends
		friend class math::SpanX;
		friend class math::MatrixX;
	};
}