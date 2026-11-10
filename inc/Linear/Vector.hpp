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
	class Span;
	class Vec3;
	class Matrix;
}

namespace math
{
	class Vector
	{
	public:
		//constructors
		Vector(const Vec3&);
		Vector(const Vector&);
		explicit Vector(uint32_t);
		Vector(double*, uint32_t);
		Vector(const double*, uint32_t);
		Vector(std::initializer_list<double>);

		//destructor
		virtual ~Vector(void);

		//size
		uint32_t rows(void) const;

		//data
		double* data(void);
		const double* data(void) const;

		//direction
		double norm(void) const;
		Vector& normalize(void);
		Vector unit(void) const;
		
		//utility
		Vector& zeros(void);
		Vector& randu(double = - 1, double = +1);

		//print
		void print(const char* = "", double = 0) const;

		//products
		Matrix outer(void) const;
		Matrix outer(const Vector&) const;
		double inner(const Vector&) const;
		double inner(const double*) const;

		//operators
		Vector operator+(void) const;
		Vector operator-(void) const;
		Vector operator/(double) const;

		Vector operator+(const Vector&) const;
		Vector operator-(const Vector&) const;

		Vector& operator=(double);
		Vector& operator=(const double*);
		Vector& operator=(const Vector&);
		Vector& operator=(std::initializer_list<double>);

		Vector& operator+=(double);
		Vector& operator-=(double);
		Vector& operator*=(double);
		Vector& operator/=(double);

		Vector& operator+=(const double*);
		Vector& operator-=(const double*);
		Vector& operator+=(const Vector&);
		Vector& operator-=(const Vector&);

		double& operator[](uint32_t);
		const double& operator[](uint32_t) const;

		friend Vector operator*(double, const Vector&);

	private:
		//data
		bool m_own;
		uint32_t m_rows;
		double* m_data_ptr;
		const double* m_data_ref;
		double m_data_mem[VECTOR_STATIC_SIZE];

		//friends
		friend class math::Span;
		friend class math::Vec3;
		friend class math::Matrix;
	};
}