//std
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>

//math
#include "Math/inc/linear/quat.hpp"
#include "Math/inc/linear/vec3.hpp"
#include "Math/inc/linear/mat3.hpp"

namespace math
{
	//constructors
	mat3::mat3(mode init) : matrix(3, 3, init)
	{
		return;
	}
	mat3::mat3(double* ptr) : matrix(ptr, 3, 3)
	{
		return;
	}
	mat3::mat3(const mat3& m) : matrix(m)
	{
		return;
	}
	mat3::mat3(const double* ref) : matrix(ref, 3, 3)
	{
		return;
	}
	mat3::mat3(std::initializer_list<double> list) : matrix(3, 3)
	{
		if(list.size() != 9)
		{
			fprintf(stderr, "Error: mat3 constructor with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		memcpy(m_ptr, std::data(list), list.size() * sizeof(double));
	}
	mat3::mat3(const vec3& v1, const vec3& v2, const vec3& v3) : matrix(3, 3)
	{
		memcpy(m_ptr + 0, v1.data(), 3 * sizeof(double));
		memcpy(m_ptr + 3, v2.data(), 3 * sizeof(double));
		memcpy(m_ptr + 6, v3.data(), 3 * sizeof(double));
	}
	mat3::mat3(const double* s1, const double* s2, const double* s3) : matrix(3, 3)
	{
		memcpy(m_ptr + 0, s1, 3 * sizeof(double));
		memcpy(m_ptr + 3, s2, 3 * sizeof(double));
		memcpy(m_ptr + 6, s3, 3 * sizeof(double));
	}

	//destructor
	mat3::~mat3(void)
	{
		return;
	}

	//operators
	mat3 mat3::operator+(void) const
	{
		return *this;
	}
	mat3 mat3::operator-(void) const
	{
		return mat3(*this) *= -1;
	}

	mat3 mat3::operator/(double s) const
	{
		return mat3(*this) /= s;
	}
	mat3 mat3::operator+(const mat3& m) const
	{
		return mat3(*this) += m;
	}
	mat3 mat3::operator-(const mat3& m) const
	{
		return mat3(*this) -= m;
	}
	mat3 mat3::operator*(const mat3& m) const
	{
		mat3 r;
		r.m_ptr[0] = m_ref[0] * m.m_ref[0] + m_ref[3] * m.m_ref[1] + m_ref[6] * m.m_ref[2];
		r.m_ptr[1] = m_ref[1] * m.m_ref[0] + m_ref[4] * m.m_ref[1] + m_ref[7] * m.m_ref[2];
		r.m_ptr[2] = m_ref[2] * m.m_ref[0] + m_ref[5] * m.m_ref[1] + m_ref[8] * m.m_ref[2];
		r.m_ptr[3] = m_ref[0] * m.m_ref[3] + m_ref[3] * m.m_ref[4] + m_ref[6] * m.m_ref[5];
		r.m_ptr[4] = m_ref[1] * m.m_ref[3] + m_ref[4] * m.m_ref[4] + m_ref[7] * m.m_ref[5];
		r.m_ptr[5] = m_ref[2] * m.m_ref[3] + m_ref[5] * m.m_ref[4] + m_ref[8] * m.m_ref[5];
		r.m_ptr[6] = m_ref[0] * m.m_ref[6] + m_ref[3] * m.m_ref[7] + m_ref[6] * m.m_ref[8];
		r.m_ptr[7] = m_ref[1] * m.m_ref[6] + m_ref[4] * m.m_ref[7] + m_ref[7] * m.m_ref[8];
		r.m_ptr[8] = m_ref[2] * m.m_ref[6] + m_ref[5] * m.m_ref[7] + m_ref[8] * m.m_ref[8];
		return r;
	}
	vec3 mat3::operator*(const vec3& v) const
	{
		vec3 r;
		r[0] = m_ref[0] * v[0] + m_ref[3] * v[1] + m_ref[6] * v[2];
		r[1] = m_ref[1] * v[0] + m_ref[4] * v[1] + m_ref[7] * v[2];
		r[2] = m_ref[2] * v[0] + m_ref[5] * v[1] + m_ref[8] * v[2];
		return r;
	}

	mat3& mat3::operator+=(double s)
	{
		for(uint32_t i = 0; i < 9; i++)
		{
			m_ptr[i] += s;
		}
		return *this;
	}
	mat3& mat3::operator-=(double s)
	{
		for(uint32_t i = 0; i < 9; i++)
		{
			m_ptr[i] -= s;
		}
		return *this;
	}
	mat3& mat3::operator*=(double s)
	{
		for(uint32_t i = 0; i < 9; i++)
		{
			m_ptr[i] *= s;
		}
		return *this;
	}
	mat3& mat3::operator/=(double s)
	{
		for(uint32_t i = 0; i < 9; i++)
		{
			m_ptr[i] /= s;
		}
		return *this;
	}

	mat3& mat3::operator=(const mat3& m)
	{
		memcpy(m_ptr, m.m_ref, 9 * sizeof(double));
		return *this;
	}

	mat3& mat3::operator+=(const mat3& m)
	{
		for(uint32_t i = 0; i < 9; i++)
		{
			m_ptr[i] += m.m_ref[i];
		}
		return *this;
	}
	mat3& mat3::operator-=(const mat3& m)
	{
		for(uint32_t i = 0; i < 9; i++)
		{
			m_ptr[i] -= m.m_ref[i];
		}
		return *this;
	}

	double& mat3::operator[](uint32_t i)
	{
		return m_ptr[i];
	}
	double& mat3::operator()(uint32_t i)
	{
		return m_ptr[i];
	}
	double& mat3::operator()(uint32_t i, uint32_t j)
	{
		return m_ptr[i + 3 * j];
	}

	const double& mat3::operator[](uint32_t i) const
	{
		return m_ref[i];
	}
	const double& mat3::operator()(uint32_t i) const
	{
		return m_ref[i];
	}
	const double& mat3::operator()(uint32_t i, uint32_t j) const
	{
		return m_ref[i + 3 * j];
	}

	//linear
	mat3 mat3::eye(void)
	{
		mat3 m;
		m.matrix::eye();
		return m;
	}
	vec3 mat3::eigen(void) const
	{
		vec3 s;
		const double tl = lode_angle();
		const double I1 = invariant_I1();
		const double J2 = deviatoric_J2();
		s[2] = I1 / 3 + 2 * sqrt(J2 / 3) * cos(tl);
		s[0] = I1 / 3 + 2 * sqrt(J2 / 3) * cos(tl + 2 * M_PI / 3);
		s[1] = I1 / 3 + 2 * sqrt(J2 / 3) * cos(tl - 2 * M_PI / 3);
		return s;
	}
	mat3 mat3::inverse(void) const
	{
		mat3 r;
		const double d = determinant();
		r.m_ptr[0] = (m_ref[4] * m_ref[8] - m_ref[5] * m_ref[7]) / d;
		r.m_ptr[1] = (m_ref[2] * m_ref[7] - m_ref[1] * m_ref[8]) / d;
		r.m_ptr[2] = (m_ref[1] * m_ref[5] - m_ref[2] * m_ref[4]) / d;
		r.m_ptr[3] = (m_ref[5] * m_ref[6] - m_ref[3] * m_ref[8]) / d;
		r.m_ptr[4] = (m_ref[0] * m_ref[8] - m_ref[2] * m_ref[6]) / d;
		r.m_ptr[5] = (m_ref[2] * m_ref[3] - m_ref[0] * m_ref[5]) / d;
		r.m_ptr[6] = (m_ref[3] * m_ref[7] - m_ref[4] * m_ref[6]) / d;
		r.m_ptr[7] = (m_ref[6] * m_ref[1] - m_ref[0] * m_ref[7]) / d;
		r.m_ptr[8] = (m_ref[0] * m_ref[4] - m_ref[1] * m_ref[3]) / d;
		return r;
	}
	mat3 mat3::transpose(void) const
	{
		mat3 r;
		r.m_ptr[0] = m_ref[0];
		r.m_ptr[1] = m_ref[3];
		r.m_ptr[2] = m_ref[6];
		r.m_ptr[3] = m_ref[1];
		r.m_ptr[4] = m_ref[4];
		r.m_ptr[5] = m_ref[7];
		r.m_ptr[6] = m_ref[2];
		r.m_ptr[7] = m_ref[5];
		r.m_ptr[8] = m_ref[8];
		return r;
	}

	double mat3::trace(void) const
	{
		return m_ref[0] + m_ref[4] + m_ref[8];
	}
	double mat3::lode_angle(void) const
	{
		const double J2 = deviatoric_J2();
		const double J3 = deviatoric_J3();
		return fabs(J2) > 1e-5 ? acos(3 * sqrt(3) / 2 * J3 / sqrt(J2 * J2 * J2)) / 3 : 0;
	}
	double mat3::determinant(void) const
	{
		return
			m_ref[0] * (m_ref[4] * m_ref[8] - m_ref[5] * m_ref[7]) +
			m_ref[3] * (m_ref[2] * m_ref[7] - m_ref[1] * m_ref[8]) +
			m_ref[6] * (m_ref[1] * m_ref[5] - m_ref[2] * m_ref[4]);
	}
	double mat3::invariant_I1(void) const
	{
		return m_ref[0] + m_ref[4] + m_ref[8];
	}
	double mat3::invariant_I2(void) const
	{
		return
			m_ref[0] * m_ref[4] + m_ref[0] * m_ref[8] + m_ref[4] * m_ref[8] -
			m_ref[1] * m_ref[1] - m_ref[2] * m_ref[2] - m_ref[5] * m_ref[5];
	}
	double mat3::invariant_I3(void) const
	{
		return
			m_ref[0] * m_ref[4] * m_ref[8] - 2 * m_ref[1] * m_ref[2] * m_ref[5] -
			m_ref[0] * m_ref[5] * m_ref[5] - m_ref[4] * m_ref[2] * m_ref[2] - m_ref[8] * m_ref[1] * m_ref[1];
	}
	double mat3::deviatoric_J2(void) const
	{
		const double I1 = invariant_I1();
		const double I2 = invariant_I2();
		return I1 * I1 / 3 - I2;
	}
	double mat3::deviatoric_J3(void) const
	{
		const double I1 = invariant_I1();
		const double I2 = invariant_I2();
		const double I3 = invariant_I3();
		return 2 * I1 * I1 * I1 / 27 - I1 * I2 / 3 + I3;
	}

	//rotation
	vec3 mat3::rotation(void) const
	{
		return quaternion().pseudo();
	}
	quat mat3::quaternion(void) const
	{
		quat q;
		const double d0 = m_ref[0];
		const double d1 = m_ref[4];
		const double d2 = m_ref[8];
		const double tr = d0 + d1 + d2;
		if(tr >= d0 && tr >= d1 && tr >= d2)
		{
			q[0] = sqrt(1 + tr) / 2;
			q[1] = (m_ref[2 + 3 * 1] - m_ref[1 + 3 * 2]) / (4 * q[0]);
			q[2] = (m_ref[0 + 3 * 2] - m_ref[2 + 3 * 0]) / (4 * q[0]);
			q[3] = (m_ref[1 + 3 * 0] - m_ref[0 + 3 * 1]) / (4 * q[0]);
		}
		else if(d0 >= d1 && d0 >= d2)
		{
			q[1] = sqrt(1 - tr + 2 * d0) / 2;
			q[0] = (m_ref[2 + 3 * 1] - m_ref[1 + 3 * 2]) / (4 * q[1]);
			q[2] = (m_ref[0 + 3 * 1] + m_ref[1 + 3 * 0]) / (4 * q[1]);
			q[3] = (m_ref[0 + 3 * 2] + m_ref[2 + 3 * 0]) / (4 * q[1]);
		}
		else if(d1 >= d2)
		{
			q[2] = sqrt(1 - tr + 2 * d1) / 2;
			q[0] = (m_ref[0 + 3 * 2] - m_ref[2 + 3 * 0]) / (4 * q[2]);
			q[1] = (m_ref[0 + 3 * 1] + m_ref[1 + 3 * 0]) / (4 * q[2]);
			q[3] = (m_ref[1 + 3 * 2] + m_ref[2 + 3 * 1]) / (4 * q[2]);
		}
		else
		{
			q[3] = sqrt(1 - tr + 2 * d2) / 2;
			q[0] = (m_ref[1 + 3 * 0] - m_ref[0 + 3 * 1]) / (4 * q[3]);
			q[1] = (m_ref[0 + 3 * 2] + m_ref[2 + 3 * 0]) / (4 * q[3]);
			q[2] = (m_ref[1 + 3 * 2] + m_ref[2 + 3 * 1]) / (4 * q[3]);
		}
		return q;
	}

	//friends
	mat3 operator*(double s, const mat3& m)
	{
		return mat3(m) *= s;
	}
}