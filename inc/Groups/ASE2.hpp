#pragma once

//Math
#include "Math/inc/Groups/LieAlgebra.hpp"

namespace math
{
	class GSE2;
}

namespace math
{
	class ASE2 : public LieAlgebra<GSE2, ASE2>
	{
	public:
		//constructor
		explicit ASE2(const double*);
		explicit ASE2(double, double, double);

		//destructor
		~ASE2(void) = default;

		//group
		GSE2 exp(void) const override;
		ASE2 bracket(const ASE2&) const override;

		//linear
		Vector vector(void) const override;
		Matrix matrix(void) const override;
		Matrix adjoint(void) const override;
		Matrix jacobian(void) const override;
		Matrix jacobian_inverse(void) const override;

		//operators
		ASE2& operator*=(double) override;
		ASE2& operator/=(double) override;
		ASE2& operator+=(const ASE2&) override;
		ASE2& operator-=(const ASE2&) override;

		ASE2 operator+(void) const override;
		ASE2 operator-(void) const override;
		ASE2 operator*(double) const override;
		ASE2 operator/(double) const override;
		ASE2 operator+(const ASE2&) const override;
		ASE2 operator-(const ASE2&) const override;

		friend ASE2 operator*(double, const ASE2&);

	private:
		//data
		double m_data[3];
	};
}