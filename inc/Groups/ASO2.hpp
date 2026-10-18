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
		explicit ASO2(const double*);

		//destructor
		~ASO2(void) = default;

		//group
		GSO2 exp(void) const override;
		ASO2 bracket(const ASO2&) const override;

		//random
		ASO2& randu(double, double) override;

		//linear
		Vector vector(void) override;
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

		friend ASO2 operator*(double, const ASO2&);

		double& operator[](uint32_t) override;
		const double& operator[](uint32_t) const override;

		//tests
		static bool test_adjoint(uint32_t, double);
		static bool test_bracket(uint32_t, double);
		static bool test_exp_log(uint32_t, double);
		static bool test_jacobian(uint32_t, double);
		static bool test_jacobian_inverse(uint32_t, double);

	private:
		//data
		double m_data;
	};
}