#pragma once

//Math
#include "Math/inc/Groups/LieAlgebra.hpp"

namespace math
{
	class GSE3;
}

namespace math
{
	class ASE3 : public LieAlgebra<GSE3, ASE3>
	{
	public:
		//constructor
		explicit ASE3(const double*);
		explicit ASE3(double = 0, double = 0, double = 0, double = 0, double = 0, double = 0);

		//destructor
		~ASE3(void) = default;

		//group
		GSE3 exp(void) const override;
		ASE3 bracket(const ASE3&) const override;

		//random
		ASE3& randu(double, double) override;

		//linear
		Vector vector(void) override;
		Vector vector(void) const override;
		Matrix matrix(void) const override;
		Matrix adjoint(void) const override;
		Matrix jacobian(void) const override;
		Matrix jacobian_inverse(void) const override;

		//operators
		ASE3& operator*=(double) override;
		ASE3& operator/=(double) override;
		ASE3& operator+=(const ASE3&) override;
		ASE3& operator-=(const ASE3&) override;

		ASE3 operator+(void) const override;
		ASE3 operator-(void) const override;
		ASE3 operator*(double) const override;
		ASE3 operator/(double) const override;
		ASE3 operator+(const ASE3&) const override;
		ASE3 operator-(const ASE3&) const override;

		friend ASE3 operator*(double, const ASE3&);

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
		double m_data[6];
	};
}