//std
#include <cmath>

//Math
#include "Math/inc/Groups/ASO2.hpp"
#include "Math/inc/Groups/GSO2.hpp"
#include "Math/inc/Miscellaneous/util.hpp"

namespace math
{
	//constructor
	GSO2::GSO2(double data) : m_data{data}
	{
		return;
	}

	//random
	GSO2& GSO2::randu(double a, double b)
	{
		return *this = ASO2().randu(a, b).exp();
	}

	//actions
	Matrix GSO2::matrix(void) const
	{
		const double c = cos(m_data);
		const double s = sin(m_data);
		return Matrix({{c, -s}, {s, c}});
	}
	Vector GSO2::action(const Vector& v) const
	{
		const double c = cos(m_data);
		const double s = sin(m_data);
		return Vector({c * v[0] - s * v[1], s * v[0] + c * v[1]});
	}

	//operations
	ASO2 GSO2::log(void) const
	{
		return ASO2(m_data);
	}
	GSO2 GSO2::inverse(void) const
	{
		return GSO2(-m_data);
	}
	GSO2 GSO2::operator*(const GSO2& element) const
	{
		return GSO2(m_data + element.m_data);
	}
	GSO2 GSO2::interpolate(const GSO2& element, double t) const
	{
		return GSO2((1 - t) * m_data + t * element.m_data);
	}

	//tests
	bool GSO2::test_matrix(uint32_t nt, double tolerance)
	{
		//test
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			const GSO2 g1(math::randu(-1, +1));
			const GSO2 g2(math::randu(-1, +1));
			test = test && ((g1 * g2).matrix() - g1.matrix() * g2.matrix()).norm() < tolerance;
		}
		//print
		printf("GSO2 matrix test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool GSO2::test_inverse(uint32_t nt, double tolerance)
	{
		//test
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			const GSO2 g(math::randu(-1, +1));
			test = test && (g * g.inverse()).log().vector().norm() < tolerance;
		}
		//print
		printf("GSO2 inverse test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
	bool GSO2::test_log_exp(uint32_t nt, double tolerance)
	{
		//test
		bool test = true;
		for(uint32_t i = 0; i < nt; i++)
		{
			const GSO2 g(math::randu(-1, +1));
			test = test && (g.matrix() - g.log().exp().matrix()).norm() < tolerance;
		}
		//print
		printf("GSO2 log-exp test:\t\t%s!\n", test ? "ok" : "not ok");
		//return
		return test;
	}
}