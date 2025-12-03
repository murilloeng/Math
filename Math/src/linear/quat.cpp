//std
#include <cmath>

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/linear/vec3.hpp"
#include "Math/Math/inc/linear/mat3.hpp"
#include "Math/Math/inc/linear/quat.hpp"

namespace math
{
	//constructors
	quat::quat(void) : vector(4)
	{
		return;
	}
	quat::quat(double* ptr) : vector(ptr, 4)
	{
		return;
	}
	quat::quat(const quat& q) : vector(q)
	{
		return;
	}
	quat::quat(const double* ref) : vector(ref, 4)
	{
		return;
	}
	quat::quat(double t, uint32_t i) : vector(4)
	{
		m_data_ptr[0] = cos(t / 2);
		m_data_ptr[1] = i == 0 ? sin(t / 2) : 0;
		m_data_ptr[2] = i == 1 ? sin(t / 2) : 0;
		m_data_ptr[3] = i == 2 ? sin(t / 2) : 0;
	}
	quat::quat(double t, const vec3& v) : vector(4)
	{
		m_data_ptr[0] = cos(t / 2);
		m_data_ptr[1] = sin(t / 2) * v[0];
		m_data_ptr[2] = sin(t / 2) * v[1];
		m_data_ptr[3] = sin(t / 2) * v[2];
	}
	quat::quat(double v0, double v1, double v2, double v3) : vector({v0, v1, v2, v3})
	{
		return;
	}
	quat::quat(const vec3& v1, const vec3& v2, const vec3& v3) : vector(4)
	{
		*this = mat3(v1, v2, v3).quaternion();
	}
	quat::quat(const double* s1, const double* s2, const double* s3) : vector(4)
	{
		*this = mat3(s1, s2, s3).quaternion();
	}

	//destructor
	quat::~quat(void)
	{
		return;
	}

	//operators
	quat& quat::operator=(const quat& q)
	{
		m_data_ptr[0] = q.m_data_ref[0];
		m_data_ptr[1] = q.m_data_ref[1];
		m_data_ptr[2] = q.m_data_ref[2];
		m_data_ptr[3] = q.m_data_ref[3];
		return *this;
	}
	quat& quat::operator*=(const quat& q)
	{
		return *this = q * *this;
	}

	quat quat::operator*(const quat& q) const
	{
		quat r;
		r.m_data_ptr[0] = m_data_ref[0] * q.m_data_ref[0] - m_data_ref[1] * q.m_data_ref[1] - m_data_ref[2] * q.m_data_ref[2] - m_data_ref[3] * q.m_data_ref[3];
		r.m_data_ptr[1] = m_data_ref[0] * q.m_data_ref[1] + m_data_ref[1] * q.m_data_ref[0] + m_data_ref[2] * q.m_data_ref[3] - m_data_ref[3] * q.m_data_ref[2];
		r.m_data_ptr[2] = m_data_ref[0] * q.m_data_ref[2] + m_data_ref[2] * q.m_data_ref[0] + m_data_ref[3] * q.m_data_ref[1] - m_data_ref[1] * q.m_data_ref[3];
		r.m_data_ptr[3] = m_data_ref[0] * q.m_data_ref[3] + m_data_ref[3] * q.m_data_ref[0] + m_data_ref[1] * q.m_data_ref[2] - m_data_ref[2] * q.m_data_ref[1];
		return r;
	}

	double& quat::operator[](uint32_t i)
	{
		return m_data_ptr[i];
	}
	double& quat::operator()(uint32_t i)
	{
		return m_data_ptr[i];
	}

	const double& quat::operator[](uint32_t i) const
	{
		return m_data_ref[i];
	}
	const double& quat::operator()(uint32_t i) const
	{
		return m_data_ref[i];
	}

	//views
	quat& quat::reset(void)
	{
		m_data_ptr[0] = 1;
		m_data_ptr[1] = 0;
		m_data_ptr[2] = 0;
		m_data_ptr[3] = 0;
		return *this;
	}
	quat& quat::view_x(void)
	{
		return *this = quat(-M_PI / 2, 2) * quat(-M_PI / 2, 1);
	}
	quat& quat::view_y(void)
	{
		return *this = quat(+M_PI / 2, 2) * quat(+M_PI / 2, 0);
	}
	quat& quat::view_z(void)
	{
		return reset();
	}
	quat& quat::isometric(uint32_t d)
	{
		switch(d)
		{
			case 0:
				view_y();
				return *this *= quat(atan(M_SQRT1_2), 0) * quat(-M_PI / 4, 1);
			case 1:
				view_z();
				return *this *= quat(atan(M_SQRT1_2), 0) * quat(-M_PI / 4, 1);
			case 2:
				view_x();
				return *this *= quat(atan(M_SQRT1_2), 0) * quat(-M_PI / 4, 1);
			default:
				return *this;
		}
	}

	//linear
	quat& quat::normalize(void)
	{
		vector::normalize();
		return *this;
	}
	vec3 quat::axial(void) const
	{
		//data
		vec3 r;
		const double s = math::vec3(m_data_ref + 1).norm();
		//axial
		r[0] = s ? m_data_ref[1] / s : 0;
		r[1] = s ? m_data_ref[2] / s : 0;
		r[2] = s ? m_data_ref[3] / s : 1;
		return r;
	}
	vec3 quat::pseudo(void) const
	{
		return angle() * axial();
	}
	double quat::angle(void) const
	{
		return 2 * acos(bound(m_data_ref[0]));
	}
	quat& quat::randu(double a, double b)
	{
		matrix::randu(a, b);
		return normalize();
	}
	vec3 quat::pseudo(const vec3& r) const
	{
		const vec3 n = axial();
		const double t = angle();
		const double k = (n.inner(r) - t) / (2 * M_PI);
		return (t + 2 * M_PI * round(k)) * n;
	}
	vec3 quat::pseudo(const quat& q) const
	{
		return pseudo(q.pseudo());
	}

	quat quat::conjugate(void) const
	{
		quat q;
		q.m_data_ptr[0] = +m_data_ref[0];
		q.m_data_ptr[1] = -m_data_ref[1];
		q.m_data_ptr[2] = -m_data_ref[2];
		q.m_data_ptr[3] = -m_data_ref[3];
		return q;
	}

	vec3 quat::rotate(const vec3& v) const
	{
		vec3 r;
		const vec3 x(m_data_ref + 1);
		const double s = m_data_ref[0];
		const double b = 2 * x.inner(v);
		const double a = s * s - x.inner(x);
		r[0] = a * v[0] + b * x[0] + 2 * s * (x[1] * v[2] - x[2] * v[1]);
		r[1] = a * v[1] + b * x[1] + 2 * s * (x[2] * v[0] - x[0] * v[2]);
		r[2] = a * v[2] + b * x[2] + 2 * s * (x[0] * v[1] - x[1] * v[0]);
		return r;
	}
	vec3 quat::conjugate(const vec3& v) const
	{
		return conjugate().rotate(v);
	}
	quat quat::conjugate(const quat& q) const
	{
		return conjugate() * q;
	}

	quat quat::slerp(const quat& q, double s) const
	{
		const quat qr = conjugate(q);
		return *this * quat(s * qr.angle(), qr.axial());
	}

	mat3 quat::rotation(void) const
	{
		mat3 m;
		m[1 + 3 * 0] = 2 * (m_data_ref[1] * m_data_ref[2] + m_data_ref[0] * m_data_ref[3]);
		m[2 + 3 * 0] = 2 * (m_data_ref[1] * m_data_ref[3] - m_data_ref[0] * m_data_ref[2]);
		m[0 + 3 * 1] = 2 * (m_data_ref[1] * m_data_ref[2] - m_data_ref[0] * m_data_ref[3]);
		m[2 + 3 * 1] = 2 * (m_data_ref[2] * m_data_ref[3] + m_data_ref[0] * m_data_ref[1]);
		m[0 + 3 * 2] = 2 * (m_data_ref[1] * m_data_ref[3] + m_data_ref[0] * m_data_ref[2]);
		m[1 + 3 * 2] = 2 * (m_data_ref[2] * m_data_ref[3] - m_data_ref[0] * m_data_ref[1]);
		m[0 + 3 * 0] = m_data_ref[0] * m_data_ref[0] + m_data_ref[1] * m_data_ref[1] - m_data_ref[2] * m_data_ref[2] - m_data_ref[3] * m_data_ref[3];
		m[1 + 3 * 1] = m_data_ref[0] * m_data_ref[0] - m_data_ref[1] * m_data_ref[1] + m_data_ref[2] * m_data_ref[2] - m_data_ref[3] * m_data_ref[3];
		m[2 + 3 * 2] = m_data_ref[0] * m_data_ref[0] - m_data_ref[1] * m_data_ref[1] - m_data_ref[2] * m_data_ref[2] + m_data_ref[3] * m_data_ref[3];
		return m;
	}
	double* quat::rotation(double* m) const
	{
		m[3 + 4 * 3] = 1;
		m[3 + 4 * 0] = m[3 + 4 * 1] = m[3 + 4 * 2] = 0;
		m[0 + 4 * 3] = m[1 + 4 * 3] = m[2 + 4 * 3] = 0;
		m[1 + 4 * 2] = 2 * (m_data_ref[2] * m_data_ref[3] - m_data_ref[0] * m_data_ref[1]);
		m[2 + 4 * 1] = 2 * (m_data_ref[2] * m_data_ref[3] + m_data_ref[0] * m_data_ref[1]);
		m[2 + 4 * 0] = 2 * (m_data_ref[1] * m_data_ref[3] - m_data_ref[0] * m_data_ref[2]);
		m[0 + 4 * 2] = 2 * (m_data_ref[1] * m_data_ref[3] + m_data_ref[0] * m_data_ref[2]);
		m[0 + 4 * 1] = 2 * (m_data_ref[1] * m_data_ref[2] - m_data_ref[0] * m_data_ref[3]);
		m[1 + 4 * 0] = 2 * (m_data_ref[1] * m_data_ref[2] + m_data_ref[0] * m_data_ref[3]);
		m[0 + 4 * 0] = m_data_ref[0] * m_data_ref[0] + m_data_ref[1] * m_data_ref[1] - m_data_ref[2] * m_data_ref[2] - m_data_ref[3] * m_data_ref[3];
		m[1 + 4 * 1] = m_data_ref[0] * m_data_ref[0] - m_data_ref[1] * m_data_ref[1] + m_data_ref[2] * m_data_ref[2] - m_data_ref[3] * m_data_ref[3];
		m[2 + 4 * 2] = m_data_ref[0] * m_data_ref[0] - m_data_ref[1] * m_data_ref[1] - m_data_ref[2] * m_data_ref[2] + m_data_ref[3] * m_data_ref[3];
		return m;
	}
}