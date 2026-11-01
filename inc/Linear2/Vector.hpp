#pragma once

//std
#include <cstdint>
#include <initializer_list>

//Math
#include "Math/inc/Linear2/Matrix.hpp"

namespace math
{
	namespace linear
	{
		class Vector : public Matrix
		{
		public:
			//constructors
			Vector(const Vector&);
			Vector(const Matrix&);
			explicit Vector(uint32_t);
			Vector(double*, uint32_t);
			Vector(const double*, uint32_t);
			Vector(std::initializer_list<double>);
	
			//destructor
			virtual ~Vector(void) = default;
	
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
}