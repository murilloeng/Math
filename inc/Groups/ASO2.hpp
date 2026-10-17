#pragma once

//Math
#include "Math/inc/Groups/LieAlgebra.hpp"

namespace math
{
	class GSO2;
}

namespace math
{
	class ASO2 : public LieAlgebra<GSO2, ASO2>
	{
	public:
		//constructor
		explicit ASO2(double = 0);

		//destructor
		~ASO2(void) = default;

		//group
		GSO2 exp(void) const override;
		ASO2 bracket(const ASO2&) const override;

		//linear
		Vector vector(void) const override;
		Matrix matrix(void) const override;
		Matrix adjoint(void) const override;
		Matrix jacobian(void) const override;
		Matrix jacobian_inverse(void) const override;

		//operators
		ASO2& operator*=(double) override;
		ASO2& operator/=(double) override;
		ASO2& operator+=(const ASO2&) override;
		ASO2& operator-=(const ASO2&) override;

		ASO2 operator+(void) const override;
		ASO2 operator-(void) const override;
		ASO2 operator*(double) const override;
		ASO2 operator/(double) const override;
		ASO2 operator+(const ASO2&) const override;
		ASO2 operator-(const ASO2&) const override;

	private:
		//data
		double m_data;
	};
}