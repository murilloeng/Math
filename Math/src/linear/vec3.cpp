//std
#include <cmath>

//math
#include "Math/Math/inc/linear/quat.hpp"
#include "Math/Math/inc/linear/vec3.hpp"
#include "Math/Math/inc/linear/mat3.hpp"
#include "Math/Math/inc/misc/rotations.hpp"

namespace math
{
	//constructors
	vec3::vec3(void) : vector(3)
	{
		return;
	}
	vec3::vec3(double* ptr) : vector(ptr, 3)
	{
		return;
	}
	vec3::vec3(const vec3& v) : vector(v)
	{
		return;
	}
	vec3::vec3(const double* ref) : vector(ref, 3)
	{
		return;
	}
	vec3::vec3(double v0, double v1, double v2) : vector({v0, v1, v2})
	{
		return;
	}

	//destructor
	vec3::~vec3(void)
	{
		return;
	}

	//operators
	vec3 vec3::operator+(void) const
	{
		return *this;
	}
	vec3 vec3::operator-(void) const
	{
		return vec3(*this) *= -1;
	}

	vec3 vec3::operator/(double s) const
	{
		return vec3(*this) /= s;
	}
	vec3 vec3::operator+(const vec3& v) const
	{
		return vec3(*this) += v;
	}
	vec3 vec3::operator-(const vec3& v) const
	{
		return vec3(*this) -= v;
	}

	vec3& vec3::operator+=(double s)
	{
		m_data_ptr[0] += s;
		m_data_ptr[1] += s;
		m_data_ptr[2] += s;
		return *this;
	}
	vec3& vec3::operator-=(double s)
	{
		m_data_ptr[0] -= s;
		m_data_ptr[1] -= s;
		m_data_ptr[2] -= s;
		return *this;
	}
	vec3& vec3::operator*=(double s)
	{
		m_data_ptr[0] *= s;
		m_data_ptr[1] *= s;
		m_data_ptr[2] *= s;
		return *this;
	}
	vec3& vec3::operator/=(double s)
	{
		m_data_ptr[0] /= s;
		m_data_ptr[1] /= s;
		m_data_ptr[2] /= s;
		return *this;
	}

	vec3& vec3::operator=(const vec3& v)
	{
		m_data_ptr[0] = v.m_data_ref[0];
		m_data_ptr[1] = v.m_data_ref[1];
		m_data_ptr[2] = v.m_data_ref[2];
		return *this;
	}

	vec3& vec3::operator+=(const vec3& v)
	{
		m_data_ptr[0] += v.m_data_ref[0];
		m_data_ptr[1] += v.m_data_ref[1];
		m_data_ptr[2] += v.m_data_ref[2];
		return *this;
	}
	vec3& vec3::operator-=(const vec3& v)
	{
		m_data_ptr[0] -= v.m_data_ref[0];
		m_data_ptr[1] -= v.m_data_ref[1];
		m_data_ptr[2] -= v.m_data_ref[2];
		return *this;
	}
	vec3& vec3::operator*=(const mat3& m)
	{
		return *this = m * *this;
	}

	//linear
	vec3& vec3::normalize(void)
	{
		return *this /= norm();
	}
	vec3& vec3::project(const vec3& v)
	{
		return *this -= inner(v) * v;
	}
	const vec3& vec3::triad(vec3& t2, vec3& t3, double a) const
	{
		//data
		uint32_t i;
		min(true, &i);
		//axes 2
		t2[(i + 0) % 3] = 0;
		t2[(i + 1) % 3] = -m_data_ref[(i + 2) % 3];
		t2[(i + 2) % 3] = +m_data_ref[(i + 1) % 3];
		//axes 3
		t3 = cross(t2.normalize());
		//rotate
		t2 = quat(a, *this).rotate(t2);
		t3 = quat(a, *this).rotate(t3);
		//return
		return *this;
	}

	mat3 vec3::spin(void) const
	{
		mat3 s;
		s[7] = -m_data_ref[0];
		s[2] = -m_data_ref[1];
		s[3] = -m_data_ref[2];
		s[5] = +m_data_ref[0];
		s[6] = +m_data_ref[1];
		s[1] = +m_data_ref[2];
		s[0] = s[4] = s[8] = 0;
		return s;
	}
	quat vec3::quaternion(void) const
	{
		//angle
		const double t = norm();
		const double c = cos(t / 2);
		const double s = t ? sin(t / 2) / t : 0.5;
		//setup
		quat q;
		q[0] = c;
		q[1] = s * m_data_ref[0];
		q[2] = s * m_data_ref[1];
		q[3] = s * m_data_ref[2];
		return q;
	}
	mat3 vec3::projection(void) const
	{
		mat3 p;
		p[0] = 1 - m_data_ref[0] * m_data_ref[0];
		p[4] = 1 - m_data_ref[1] * m_data_ref[1];
		p[8] = 1 - m_data_ref[2] * m_data_ref[2];
		p[1] = p[3] = -m_data_ref[0] * m_data_ref[1];
		p[2] = p[6] = -m_data_ref[0] * m_data_ref[2];
		p[5] = p[7] = -m_data_ref[1] * m_data_ref[2];
		return p;
	}

	double vec3::inner(const vec3& v) const
	{
		return
			m_data_ref[0] * v.m_data_ref[0] +
			m_data_ref[1] * v.m_data_ref[1] +
			m_data_ref[2] * v.m_data_ref[2];
	}
	mat3 vec3::outer(const vec3& v) const
	{
		mat3 m;
		m[0] = m_data_ref[0] * v.m_data_ref[0];
		m[1] = m_data_ref[1] * v.m_data_ref[0];
		m[2] = m_data_ref[2] * v.m_data_ref[0];
		m[3] = m_data_ref[0] * v.m_data_ref[1];
		m[4] = m_data_ref[1] * v.m_data_ref[1];
		m[5] = m_data_ref[2] * v.m_data_ref[1];
		m[6] = m_data_ref[0] * v.m_data_ref[2];
		m[7] = m_data_ref[1] * v.m_data_ref[2];
		m[8] = m_data_ref[2] * v.m_data_ref[2];
		return m;
	}
	vec3 vec3::cross(const vec3& v) const
	{
		vec3 r;
		r.m_data_ptr[0] = m_data_ref[1] * v.m_data_ref[2] - m_data_ref[2] * v.m_data_ref[1];
		r.m_data_ptr[1] = m_data_ref[2] * v.m_data_ref[0] - m_data_ref[0] * v.m_data_ref[2];
		r.m_data_ptr[2] = m_data_ref[0] * v.m_data_ref[1] - m_data_ref[1] * v.m_data_ref[0];
		return r;
	}
	vec3 vec3::rotate(const vec3& v) const
	{
		const double t = norm();
		return v + fn(t, 1) * cross(v) + fn(t, 2) * cross(cross(v));
	}

	vec3 vec3::lerp(const vec3& v, double s) const
	{
		return (1 - s) * *this + s * v;
	}

	mat3 vec3::rotation_tensor(void) const
	{
		const double t = norm();
		return mat3::eye() + fn(t, 1) * spin() + fn(t, 2) * spin() * spin();
	}

	mat3 vec3::rotation_gradient(bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		return mat3::eye() + s * fn(t, 2) * spin() + fn(t, 3) * spin() * spin();
	}
	vec3 vec3::rotation_gradient(const vec3& u, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		return u + s * fn(t, 2) * cross(u) + fn(t, 3) * cross(cross(u));
	}

	mat3 vec3::rotation_class(uint32_t n, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		const double f0 = fn(0, n + 0);
		const double f1 = fn(t, n + 1);
		const double f2 = fn(t, n + 2);
		return f0 * mat3::eye() + s * f1 * spin() + f2 * spin() * spin();
	}
	vec3 vec3::rotation_class(const vec3& v, uint32_t n, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		const double f0 = fn(0, n + 0);
		const double f1 = fn(t, n + 1);
		const double f2 = fn(t, n + 2);
		return f0 * v + s * f1 * cross(v) + f2 * cross(cross(v));
	}

	mat3 vec3::rotation_gradient_inverse(bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		return mat3::eye() - s * spin() / 2 + (fn(t, 3) - 2 * fn(t, 4)) / fn(t, 2) / 2 * spin() * spin();
	}
	vec3 vec3::rotation_gradient_inverse(const vec3& u, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		return u - s * cross(u) / 2 + (fn(t, 3) - 2 * fn(t, 4)) / fn(t, 2) / 2 * cross(cross(u));
	}

	mat3 vec3::rotation_hessian(const vec3& u, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		const double a = 2 * fn(t, 4) - fn(t, 3);
		const double b = 3 * fn(t, 5) - fn(t, 4);
		return s * a * cross(u).outer(*this) + b * cross(cross(u)).outer(*this) - s * fn(t, 2) * u.spin() - fn(t, 3) * (cross(u).spin() + spin() * u.spin());
	}
	vec3 vec3::rotation_hessian(const vec3& u, const vec3& v, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		const double a = 2 * fn(t, 4) - fn(t, 3);
		const double b = 3 * fn(t, 5) - fn(t, 4);
		return s * a * inner(v) * cross(u) + b * inner(v) * cross(cross(u)) - s * fn(t, 2) * u.cross(v) - fn(t, 3) * (cross(u).cross(v) + cross(u.cross(v)));
	}

	mat3 vec3::rotation_hessian_inverse(const vec3& u, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		const double a = (fn(t, 3) - 2 * fn(t, 4)) / fn(t, 2) / 2;
		const double b = (fn(t, 5) - 4 * fn(t, 6)) / fn(t, 2) / 2;
		return s * u.spin() / 2 - a * (cross(u).spin() + spin() * u.spin()) + b * cross(cross(u)).outer(*this);
	}
	vec3 vec3::rotation_hessian_inverse(const vec3& u, const vec3& v, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		const double a = (fn(t, 3) - 2 * fn(t, 4)) / fn(t, 2) / 2;
		const double b = (fn(t, 5) - 4 * fn(t, 6)) / fn(t, 2) / 2;
		return s * u.cross(v) / 2 - a * (cross(u).cross(v) + cross(u.cross(v))) + b * inner(v) * cross(cross(u));
	}

	mat3 vec3::rotation_class_increment(const vec3& u, uint32_t n, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		const double f1 = fn(t, n + 1);
		const double f2 = fn(t, n + 2);
		const double df1 = (n + 1) * fn(t, n + 3) - fn(t, n + 2);
		const double df2 = (n + 2) * fn(t, n + 4) - fn(t, n + 3);
		return (s * df1 * cross(u) + df2 * cross(cross(u))).outer(*this) - s * f1 * u.spin() - f2 * (cross(u).spin() + spin() * u.spin());
	}
	vec3 vec3::rotation_class_increment(const vec3& u, const vec3& v, uint32_t n, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		const double f1 = fn(t, n + 1);
		const double f2 = fn(t, n + 2);
		const double df1 = (n + 1) * fn(t, n + 3) - fn(t, n + 2);
		const double df2 = (n + 2) * fn(t, n + 4) - fn(t, n + 3);
		return inner(v) * (s * df1 * cross(u) + df2 * cross(cross(u))) - s * f1 * u.cross(v) - f2 * (cross(u).cross(v) + cross(u.cross(v)));
	}

	mat3 vec3::rotation_third(const vec3& u, const vec3& v, bool q, bool x) const
	{
		if(x)
		{
			const double t = norm();
			const int32_t s = q ? -1 : +1;
			const double a = 2 * fn(t, 4) - fn(t, 3);
			const double b = 3 * fn(t, 5) - fn(t, 4);
			const double da = 8 * fn(t, 6) - 5 * fn(t, 5) + fn(t, 4);
			const double db = 15 * fn(t, 7) - 7 * fn(t, 6) + fn(t, 5);
			return s * da * inner(v) * cross(u).outer(*this) + s * a * (cross(u).outer(v) - inner(v) * u.spin()) +
			db * inner(v) * cross(cross(u)).outer(*this) + b * cross(cross(u)).outer(v) - b * inner(v) * (spin() * u.spin() + cross(u).spin()) -
			s * a * u.cross(v).outer(*this) - b * (cross(u).cross(v) + cross(u.cross(v))).outer(*this) - fn(t, 3) * (v.spin() * u.spin() - u.cross(v).spin());
		}
		else
		{
			const double t = norm();
			const int32_t s = q ? -1 : +1;
			const double a = 2 * fn(t, 4) - fn(t, 3);
			const double b = 3 * fn(t, 5) - fn(t, 4);
			return inner(v) * (s * a * spin() + b * spin() * spin()) + s * fn(t, 2) * v.spin() + fn(t, 3) * (v.spin() * spin() + spin() * v.spin());
		}
	}
	mat3 vec3::rotation_third_inverse(const vec3& u, const vec3& v, bool q, bool flag) const
	{
		//data
		const mat3 S = spin();
		const mat3 U = u.spin();
		const mat3 V = v.spin();
		const vec3 w = u.cross(v);
		//third
		if(flag)
		{
			const double t = norm();
			const double a = (fn(t, 3) - 2 * fn(t, 4)) / fn(t, 2) / 2;
			const double b = (fn(t, 5) - 4 * fn(t, 6)) / fn(t, 2) / 2;
			const double c = (9 * fn(t, 7) - fn(t, 6) - 24 * fn(t, 8)) / fn(t, 2) / 2 + 2 * a * b;
			return
				a * (U * V - 2 * V * U) -
				b * (cross(u).cross(v) + cross(u.cross(v))).outer(*this) +
				b * inner(v) * (U * S - 2 * S * U) + b * cross(cross(u)).outer(v) + 
				c * inner(v) * cross(cross(u)).outer(*this);
		}
		else
		{
			const mat3 S = spin();
			const double t = norm();
			const int32_t s = q ? -1 : +1;
			const double a = (fn(t, 3) - 2 * fn(t, 4)) / fn(t, 2) / 2;
			const double b = (fn(t, 5) - 4 * fn(t, 6)) / fn(t, 2) / 2;
			return -s * V / 2 + a * (V * S + S * V) + b * inner(v) * S * S;
		}
	}

	//friends
	vec3 operator*(double s, const vec3& v)
	{
		return vec3(v) *= s;
	}
}