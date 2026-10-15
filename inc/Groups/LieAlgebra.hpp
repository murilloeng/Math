#pragma once

//std
#include <cstdint>

//Math
#include "Math/inc/Linear/Matrix.hpp"

namespace math
{
	class LieGroup;
}

namespace math
{
	class LieAlgebra
	{
	public:
		//constructor
		LieAlgebra(void);

		//destructor
		~LieAlgebra(void);

		//data
		double* data(void);
		const double* data(void) const;
		uint32_t dimension(void) const;

		//group
		virtual void exp(LieGroup&) const = 0;
		virtual Matrix matrix(void) const = 0;
		virtual Matrix adjoint(void) const = 0;
		virtual Matrix gradient(void) const = 0;
		virtual LieAlgebra bracket(const LieAlgebra&) const = 0;

		//operators
		LieAlgebra& operator*=(double);
		LieAlgebra& operator/=(double);
		LieAlgebra& operator+=(const LieAlgebra&);
		LieAlgebra& operator-=(const LieAlgebra&);

		LieAlgebra operator+(void) const;
		LieAlgebra operator-(void) const;
		LieAlgebra operator*(double) const;
		LieAlgebra operator/(double) const;
		LieAlgebra operator+(const LieAlgebra&) const;
		LieAlgebra operator-(const LieAlgebra&) const;

		friend LieAlgebra operator*(double, const LieAlgebra&);

	private:
		//data
		double m_data[6];
		uint32_t m_dimension;
	};
}