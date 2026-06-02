#pragma once

//std
#include <cstdint>
#include <initializer_list>

//Math
#include "Math/inc/Linear/Matrix.hpp"

namespace math
{
	class Vector : public Matrix
	{
	public:
		//constructors
		Vector(void);
		Vector(const Matrix&);
		Vector(const double*, uint32_t);
		Vector(uint32_t, mode = mode::null);
		Vector(std::initializer_list<double>);
		Vector(double*, uint32_t, mode = mode::null);

		//destructor
		virtual ~Vector(void);

		//size
		Vector& resize(uint32_t);
		Vector& resize(uint32_t, uint32_t);

		//linear
		Vector& normalize(void);
		Vector unit(void) const;
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
	};
}