#pragma once

//Math
#include "Math/inc/Groups/LieAlgebra.hpp"

namespace math
{
	class GSO3;
}

namespace math
{
	class ASO3 : public LieAlgebra<GSO3, ASO3, 3, 3>
	{
	public:
		//constructor
		explicit ASO3(const Vec3&);
		explicit ASO3(const double*);
		explicit ASO3(double = 0, double = 0, double = 0);

		//destructor
		~ASO3(void) = default;

		//group
		GSO3 exp(void) const override;
		ASO3 bracket(const ASO3&) const override;

		//random
		ASO3& randu(double, double) override;

		//linear
		Vec3 vector(void) const override;
		Mat3 matrix(void) const override;
		Mat3 adjoint(void) const override;
		Mat3 jacobian(void) const override;
		Mat3 jacobian_inverse(void) const override;
		Mat3 jacobian_gradient(const Vec3&) const override;
		Mat3 jacobian_inverse_gradient(const Vec3&, bool) const;
		Mat3 jacobian_inverse_hessian(const Vec3&, const Vec3&, bool, bool) const;

		//operators
		ASO3& operator*=(double) override;
		ASO3& operator/=(double) override;
		ASO3& operator+=(const ASO3&) override;
		ASO3& operator-=(const ASO3&) override;

		ASO3 operator+(void) const override;
		ASO3 operator-(void) const override;
		ASO3 operator*(double) const override;
		ASO3 operator/(double) const override;
		ASO3 operator+(const ASO3&) const override;
		ASO3 operator-(const ASO3&) const override;

		friend ASO3 operator*(double, const ASO3&);

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
		double m_data[3];
	};
}