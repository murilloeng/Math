//std
#include <cmath>

//Math
#include "Math/inc/Miscellaneous/util.hpp"
#include "Math/inc/Linear/Vec3.hpp"
#include "Math/inc/Linear/Mat3.hpp"
#include "Math/inc/Linear/Quat.hpp"

namespace math
{
	//constructors
	Quat::Quat(void) : Vector(4)
	{
		return;
	}
	Quat::Quat(double* ptr) : Vector(ptr, 4)
	{
		return;
	}
	Quat::Quat(const Quat& q) : Vector(q)
	{
		return;
	}
	Quat::Quat(const double* ref) : Vector(ref, 4)
	{
		return;
	}
	Quat::Quat(double t, uint32_t i) : Vector(4)
	{
		m_data_ptr[0] = cos(t / 2);
		m_data_ptr[1] = i == 0 ? sin(t / 2) : 0;
		m_data_ptr[2] = i == 1 ? sin(t / 2) : 0;
		m_data_ptr[3] = i == 2 ? sin(t / 2) : 0;
	}
	Quat::Quat(double t, const Vec3& v) : Vector(4)
	{
		m_data_ptr[0] = cos(t / 2);
		m_data_ptr[1] = sin(t / 2) * v[0];
		m_data_ptr[2] = sin(t / 2) * v[1];
		m_data_ptr[3] = sin(t / 2) * v[2];
	}
	Quat::Quat(double v0, double v1, double v2, double v3) : Vector({v0, v1, v2, v3})
	{
		return;
	}
	Quat::Quat(const Vec3& v1, const Vec3& v2, const Vec3& v3) : Vector(4)
	{
		*this = Mat3(v1, v2, v3).quaternion();
	}
	Quat::Quat(const double* s1, const double* s2, const double* s3) : Vector(4)
	{
		*this = Mat3(s1, s2, s3).quaternion();
	}

	//destructor
	Quat::~Quat(void)
	{
		return;
	}

	//operators
	Quat& Quat::operator=(const Quat& q)
	{
		m_data_ptr[0] = q.m_data_ref[0];
		m_data_ptr[1] = q.m_data_ref[1];
		m_data_ptr[2] = q.m_data_ref[2];
		m_data_ptr[3] = q.m_data_ref[3];
		return *this;
	}
	Quat& Quat::operator*=(const Quat& q)
	{
		return *this = q * *this;
	}

	Quat Quat::operator*(const Quat& q) const
	{
		Quat r;
		r.m_data_ptr[0] = m_data_ref[0] * q.m_data_ref[0] - m_data_ref[1] * q.m_data_ref[1] - m_data_ref[2] * q.m_data_ref[2] - m_data_ref[3] * q.m_data_ref[3];
		r.m_data_ptr[1] = m_data_ref[0] * q.m_data_ref[1] + m_data_ref[1] * q.m_data_ref[0] + m_data_ref[2] * q.m_data_ref[3] - m_data_ref[3] * q.m_data_ref[2];
		r.m_data_ptr[2] = m_data_ref[0] * q.m_data_ref[2] + m_data_ref[2] * q.m_data_ref[0] + m_data_ref[3] * q.m_data_ref[1] - m_data_ref[1] * q.m_data_ref[3];
		r.m_data_ptr[3] = m_data_ref[0] * q.m_data_ref[3] + m_data_ref[3] * q.m_data_ref[0] + m_data_ref[1] * q.m_data_ref[2] - m_data_ref[2] * q.m_data_ref[1];
		return r;
	}

	double& Quat::operator[](uint32_t i)
	{
		return m_data_ptr[i];
	}
	double& Quat::operator()(uint32_t i)
	{
		return m_data_ptr[i];
	}

	const double& Quat::operator[](uint32_t i) const
	{
		return m_data_ref[i];
	}
	const double& Quat::operator()(uint32_t i) const
	{
		return m_data_ref[i];
	}

	//views
	Quat& Quat::reset(void)
	{
		m_data_ptr[0] = 1;
		m_data_ptr[1] = 0;
		m_data_ptr[2] = 0;
		m_data_ptr[3] = 0;
		return *this;
	}
	Quat& Quat::view_x(void)
	{
		return *this = Quat(-M_PI / 2, 2) * Quat(-M_PI / 2, 1);
	}
	Quat& Quat::view_y(void)
	{
		return *this = Quat(+M_PI / 2, 2) * Quat(+M_PI / 2, 0);
	}
	Quat& Quat::view_z(void)
	{
		return reset();
	}
	Quat& Quat::isometric(uint32_t d)
	{
		switch(d)
		{
			case 0:
				view_y();
				return *this *= Quat(atan(M_SQRT1_2), 0) * Quat(-M_PI / 4, 1);
			case 1:
				view_z();
				return *this *= Quat(atan(M_SQRT1_2), 0) * Quat(-M_PI / 4, 1);
			case 2:
				view_x();
				return *this *= Quat(atan(M_SQRT1_2), 0) * Quat(-M_PI / 4, 1);
			default:
				return *this;
		}
	}

	//linear
	Quat& Quat::normalize(void)
	{
		Vector::normalize();
		return *this;
	}
	Vec3 Quat::axial(void) const
	{
		//data
		Vec3 r;
		const double s = math::Vec3(m_data_ref + 1).norm();
		//axial
		r[0] = s ? m_data_ref[1] / s : 0;
		r[1] = s ? m_data_ref[2] / s : 0;
		r[2] = s ? m_data_ref[3] / s : 1;
		return r;
	}
	Vec3 Quat::pseudo(void) const
	{
		//data
		const double a = m_data_ref[0];
		const double b = Vec3(m_data_ref + 1).norm();
		//check
		if(b < 1e-12) return 2 * Vec3(m_data_ref + 1);
		//return
		return 2 * atan2(b, a) / b * Vec3(m_data_ref + 1);
	}
	double Quat::angle(void) const
	{
		return 2 * acos(bound(m_data_ref[0]));
	}
	Quat& Quat::randu(double a, double b)
	{
		Matrix::randu(a, b);
		return normalize();
	}
	Vec3 Quat::pseudo(const Vec3& r) const
	{
		const Vec3 n = axial();
		const double t = angle();
		const double k = (n.inner(r) - t) / (2 * M_PI);
		return (t + 2 * M_PI * round(k)) * n;
	}
	Vec3 Quat::pseudo(const Quat& q) const
	{
		return pseudo(q.pseudo());
	}

	Quat Quat::conjugate(void) const
	{
		Quat q;
		q.m_data_ptr[0] = +m_data_ref[0];
		q.m_data_ptr[1] = -m_data_ref[1];
		q.m_data_ptr[2] = -m_data_ref[2];
		q.m_data_ptr[3] = -m_data_ref[3];
		return q;
	}

	Vec3 Quat::rotate(const Vec3& v) const
	{
		Vec3 r;
		const Vec3 x(m_data_ref + 1);
		const double s = m_data_ref[0];
		const double b = 2 * x.inner(v);
		const double a = s * s - x.inner(x);
		r[0] = a * v[0] + b * x[0] + 2 * s * (x[1] * v[2] - x[2] * v[1]);
		r[1] = a * v[1] + b * x[1] + 2 * s * (x[2] * v[0] - x[0] * v[2]);
		r[2] = a * v[2] + b * x[2] + 2 * s * (x[0] * v[1] - x[1] * v[0]);
		return r;
	}
	Vec3 Quat::conjugate(const Vec3& v) const
	{
		return conjugate().rotate(v);
	}
	Quat Quat::conjugate(const Quat& q) const
	{
		return conjugate() * q;
	}

	Quat Quat::slerp(const Quat& q, double s) const
	{
		const Quat qr = conjugate(q);
		return *this * Quat(s * qr.angle(), qr.axial());
	}

	Mat3 Quat::rotation(void) const
	{
		Mat3 m;
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
	double* Quat::rotation(double* m) const
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