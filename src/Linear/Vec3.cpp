//std
#include <cmath>

//Math
#include "Math/inc/Linear/Quat.hpp"
#include "Math/inc/Linear/Vec3.hpp"
#include "Math/inc/Linear/Mat3.hpp"
#include "Math/inc/Miscellaneous/rotations.hpp"

namespace math
{
	//constructors
	Vec3::Vec3(void) : Vector(3)
	{
		return;
	}
	Vec3::Vec3(double* ptr) : Vector(ptr, 3)
	{
		return;
	}
	Vec3::Vec3(const Vec3& v) : Vector(v)
	{
		return;
	}
	Vec3::Vec3(const double* ref) : Vector(ref, 3)
	{
		return;
	}
	Vec3::Vec3(double v0, double v1, double v2) : Vector({v0, v1, v2})
	{
		return;
	}

	//destructor
	Vec3::~Vec3(void)
	{
		return;
	}

	//operators
	Vec3 Vec3::operator+(void) const
	{
		return *this;
	}
	Vec3 Vec3::operator-(void) const
	{
		return Vec3(*this) *= -1;
	}

	Vec3 Vec3::operator/(double s) const
	{
		return Vec3(*this) /= s;
	}
	Vec3 Vec3::operator+(const Vec3& v) const
	{
		return Vec3(*this) += v;
	}
	Vec3 Vec3::operator-(const Vec3& v) const
	{
		return Vec3(*this) -= v;
	}

	Vec3& Vec3::operator+=(double s)
	{
		m_data_ptr[0] += s;
		m_data_ptr[1] += s;
		m_data_ptr[2] += s;
		return *this;
	}
	Vec3& Vec3::operator-=(double s)
	{
		m_data_ptr[0] -= s;
		m_data_ptr[1] -= s;
		m_data_ptr[2] -= s;
		return *this;
	}
	Vec3& Vec3::operator*=(double s)
	{
		m_data_ptr[0] *= s;
		m_data_ptr[1] *= s;
		m_data_ptr[2] *= s;
		return *this;
	}
	Vec3& Vec3::operator/=(double s)
	{
		m_data_ptr[0] /= s;
		m_data_ptr[1] /= s;
		m_data_ptr[2] /= s;
		return *this;
	}

	Vec3& Vec3::operator=(const Vec3& v)
	{
		m_data_ptr[0] = v.m_data_ref[0];
		m_data_ptr[1] = v.m_data_ref[1];
		m_data_ptr[2] = v.m_data_ref[2];
		return *this;
	}

	Vec3& Vec3::operator+=(const Vec3& v)
	{
		m_data_ptr[0] += v.m_data_ref[0];
		m_data_ptr[1] += v.m_data_ref[1];
		m_data_ptr[2] += v.m_data_ref[2];
		return *this;
	}
	Vec3& Vec3::operator-=(const Vec3& v)
	{
		m_data_ptr[0] -= v.m_data_ref[0];
		m_data_ptr[1] -= v.m_data_ref[1];
		m_data_ptr[2] -= v.m_data_ref[2];
		return *this;
	}
	Vec3& Vec3::operator*=(const Mat3& m)
	{
		return *this = m * *this;
	}

	//linear
	Vec3 Vec3::unit(void) const
	{
		return *this / norm();
	}
	Vec3& Vec3::normalize(void)
	{
		return *this /= norm();
	}
	Vec3& Vec3::project(const Vec3& v)
	{
		return *this -= inner(v) * v;
	}
	const Vec3& Vec3::triad(Vec3& t2, Vec3& t3, double a) const
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
		t2 = Quat(a, *this).rotate(t2);
		t3 = Quat(a, *this).rotate(t3);
		//return
		return *this;
	}

	Mat3 Vec3::spin(void) const
	{
		Mat3 s;
		s[7] = -m_data_ref[0];
		s[2] = -m_data_ref[1];
		s[3] = -m_data_ref[2];
		s[5] = +m_data_ref[0];
		s[6] = +m_data_ref[1];
		s[1] = +m_data_ref[2];
		s[0] = s[4] = s[8] = 0;
		return s;
	}
	Quat Vec3::quaternion(void) const
	{
		//angle
		const double t = norm();
		const double c = cos(t / 2);
		const double s = t ? sin(t / 2) / t : 0.5;
		//setup
		Quat q;
		q[0] = c;
		q[1] = s * m_data_ref[0];
		q[2] = s * m_data_ref[1];
		q[3] = s * m_data_ref[2];
		return q;
	}
	Mat3 Vec3::projection(void) const
	{
		Mat3 p;
		p[0] = 1 - m_data_ref[0] * m_data_ref[0];
		p[4] = 1 - m_data_ref[1] * m_data_ref[1];
		p[8] = 1 - m_data_ref[2] * m_data_ref[2];
		p[1] = p[3] = -m_data_ref[0] * m_data_ref[1];
		p[2] = p[6] = -m_data_ref[0] * m_data_ref[2];
		p[5] = p[7] = -m_data_ref[1] * m_data_ref[2];
		return p;
	}

	double Vec3::inner(const Vec3& v) const
	{
		return
			m_data_ref[0] * v.m_data_ref[0] +
			m_data_ref[1] * v.m_data_ref[1] +
			m_data_ref[2] * v.m_data_ref[2];
	}
	Mat3 Vec3::outer(const Vec3& v) const
	{
		Mat3 m;
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
	Vec3 Vec3::cross(const Vec3& v) const
	{
		Vec3 r;
		r.m_data_ptr[0] = m_data_ref[1] * v.m_data_ref[2] - m_data_ref[2] * v.m_data_ref[1];
		r.m_data_ptr[1] = m_data_ref[2] * v.m_data_ref[0] - m_data_ref[0] * v.m_data_ref[2];
		r.m_data_ptr[2] = m_data_ref[0] * v.m_data_ref[1] - m_data_ref[1] * v.m_data_ref[0];
		return r;
	}
	Vec3 Vec3::rotate(const Vec3& v) const
	{
		const double t = norm();
		return v + fn(t, 1) * cross(v) + fn(t, 2) * cross(cross(v));
	}

	Vec3 Vec3::lerp(const Vec3& v, double s) const
	{
		return (1 - s) * *this + s * v;
	}

	//static
	Vec3 Vec3::base(uint32_t index)
	{
		//data
		const double data[] = {
			1, 0, 0,
			0, 1, 0,
			0, 0, 1
		};
		const double* v = data + 3 * index;
		//return
		return Vec3(v[0], v[1], v[2]);
	}

	//rotation
	Mat3 Vec3::rotation_tensor(void) const
	{
		const double t = norm();
		return Mat3::eye() + fn(t, 1) * spin() + fn(t, 2) * spin() * spin();
	}

	Mat3 Vec3::rotation_gradient(bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		return Mat3::eye() + s * fn(t, 2) * spin() + fn(t, 3) * spin() * spin();
	}
	Vec3 Vec3::rotation_gradient(const Vec3& u, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		return u + s * fn(t, 2) * cross(u) + fn(t, 3) * cross(cross(u));
	}

	Mat3 Vec3::rotation_class(uint32_t n, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		const double f0 = fn(0, n + 0);
		const double f1 = fn(t, n + 1);
		const double f2 = fn(t, n + 2);
		return f0 * Mat3::eye() + s * f1 * spin() + f2 * spin() * spin();
	}
	Vec3 Vec3::rotation_class(const Vec3& v, uint32_t n, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		const double f0 = fn(0, n + 0);
		const double f1 = fn(t, n + 1);
		const double f2 = fn(t, n + 2);
		return f0 * v + s * f1 * cross(v) + f2 * cross(cross(v));
	}

	Mat3 Vec3::rotation_gradient_inverse(bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		return Mat3::eye() - s * spin() / 2 + (fn(t, 3) - 2 * fn(t, 4)) / fn(t, 2) / 2 * spin() * spin();
	}
	Vec3 Vec3::rotation_gradient_inverse(const Vec3& u, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		return u - s * cross(u) / 2 + (fn(t, 3) - 2 * fn(t, 4)) / fn(t, 2) / 2 * cross(cross(u));
	}

	Mat3 Vec3::rotation_hessian(uint32_t index, bool q) const
	{
		//data
		const Mat3 St = spin();
		const double t = norm();
		const double f2 = fn(t, 2);
		const double f3 = fn(t, 3);
		const int32_t s = q ? -1 : +1;
		const double tk = m_data_ref[index];
		const math::Mat3 Sk = base(index).spin();
		const double df2 = 2 * fn(t, 4) - fn(t, 3);
		const double df3 = 3 * fn(t, 5) - fn(t, 4);
		//return
		return s * df2 * tk * St + df3 * tk * St * St + s * f2 * Sk + f3 * (St * Sk + Sk * St);
	}
	Mat3 Vec3::rotation_hessian(const Vec3& u, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		const double a = 2 * fn(t, 4) - fn(t, 3);
		const double b = 3 * fn(t, 5) - fn(t, 4);
		return s * a * cross(u).outer(*this) + b * cross(cross(u)).outer(*this) - s * fn(t, 2) * u.spin() - fn(t, 3) * (cross(u).spin() + spin() * u.spin());
	}
	Vec3 Vec3::rotation_hessian(const Vec3& u, const Vec3& v, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		const double a = 2 * fn(t, 4) - fn(t, 3);
		const double b = 3 * fn(t, 5) - fn(t, 4);
		return s * a * inner(v) * cross(u) + b * inner(v) * cross(cross(u)) - s * fn(t, 2) * u.cross(v) - fn(t, 3) * (cross(u).cross(v) + cross(u.cross(v)));
	}

	Mat3 Vec3::rotation_hessian_inverse(uint32_t k, bool q) const
	{
		//data
		const double t = norm();
		const double s = q ? -1 : +1;
		const double tk = m_data_ref[k];
		const double at = (fn(t, 3) - 2 * fn(t, 4)) / fn(t, 2) / 2;
		const double bt = (fn(t, 5) - 4 * fn(t, 6)) / fn(t, 2) / 2;
		//hessian
		const Mat3 St = spin();
		const Mat3 Sk = Vec3::base(k).spin();
		return -s / 2 * Sk + at * (St * Sk + Sk * St) + bt * tk * St * St;
	}
	Mat3 Vec3::rotation_hessian_inverse(const Vec3& u, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		const double a = (fn(t, 3) - 2 * fn(t, 4)) / fn(t, 2) / 2;
		const double b = (fn(t, 5) - 4 * fn(t, 6)) / fn(t, 2) / 2;
		return s * u.spin() / 2 - a * (cross(u).spin() + spin() * u.spin()) + b * cross(cross(u)).outer(*this);
	}
	Vec3 Vec3::rotation_hessian_inverse(const Vec3& u, const Vec3& v, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		const double a = (fn(t, 3) - 2 * fn(t, 4)) / fn(t, 2) / 2;
		const double b = (fn(t, 5) - 4 * fn(t, 6)) / fn(t, 2) / 2;
		return s * u.cross(v) / 2 - a * (cross(u).cross(v) + cross(u.cross(v))) + b * inner(v) * cross(cross(u));
	}

	Mat3 Vec3::rotation_class_increment(const Vec3& u, uint32_t n, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		const double f1 = fn(t, n + 1);
		const double f2 = fn(t, n + 2);
		const double df1 = (n + 1) * fn(t, n + 3) - fn(t, n + 2);
		const double df2 = (n + 2) * fn(t, n + 4) - fn(t, n + 3);
		return (s * df1 * cross(u) + df2 * cross(cross(u))).outer(*this) - s * f1 * u.spin() - f2 * (cross(u).spin() + spin() * u.spin());
	}
	Vec3 Vec3::rotation_class_increment(const Vec3& u, const Vec3& v, uint32_t n, bool q) const
	{
		const double t = norm();
		const int32_t s = q ? -1 : +1;
		const double f1 = fn(t, n + 1);
		const double f2 = fn(t, n + 2);
		const double df1 = (n + 1) * fn(t, n + 3) - fn(t, n + 2);
		const double df2 = (n + 2) * fn(t, n + 4) - fn(t, n + 3);
		return inner(v) * (s * df1 * cross(u) + df2 * cross(cross(u))) - s * f1 * u.cross(v) - f2 * (cross(u).cross(v) + cross(u.cross(v)));
	}

	Mat3 Vec3::rotation_third(uint32_t k, uint32_t p, bool q) const
	{
		//data
		const double t = norm();
		const double f3 = fn(t, 3);
		const Vec3 ek = Vec3::base(k);
		const Vec3 ep = Vec3::base(p);
		const int32_t s = q ? -1 : +1;
		const double tk = m_data_ref[k];
		const double tp = m_data_ref[p];
		const Vec3 ea = tk * ep + tp * ek;
		const double df2 = 2 * fn(t, 4) - fn(t, 3);
		const double df3 = 3 * fn(t, 5) - fn(t, 4);
		const double ddf2 = 8 * fn(t, 6) - 5 * fn(t, 5) + fn(t, 4);
		const double ddf3 = 15 * fn(t, 7) - 7 * fn(t, 6) + fn(t, 5);
		//third
		const Mat3 St = spin();
		const Mat3 Sk = ek.spin();
		const Mat3 Sp = ep.spin();
		const Mat3 Sa = ea.spin();
		return 
			f3 * (Sk * Sp + Sp * Sk) + 
			(k == p) * (df2 * St + df3 * St * St) +
			s * df2 * Sa + df3 * (St * Sa + Sa * St) + 
			tk * tp * (s * ddf2 * St + ddf3 * St * St);
	}
	Mat3 Vec3::rotation_third(const Vec3& v, const Vec3& u, bool q, bool x) const
	{
		if(x)
		{
			const double t = norm();
			const int32_t s = q ? -1 : +1;
			const double a = 2 * fn(t, 4) - fn(t, 3);
			const double b = 3 * fn(t, 5) - fn(t, 4);
			const double da = 8 * fn(t, 6) - 5 * fn(t, 5) + fn(t, 4);
			const double db = 15 * fn(t, 7) - 7 * fn(t, 6) + fn(t, 5);
			return s * da * inner(u) * cross(v).outer(*this) + s * a * (cross(v).outer(u) - inner(u) * v.spin()) +
			db * inner(u) * cross(cross(v)).outer(*this) + b * cross(cross(v)).outer(u) - b * inner(u) * (spin() * v.spin() + cross(v).spin()) -
			s * a * v.cross(u).outer(*this) - b * (cross(v).cross(u) + cross(v.cross(u))).outer(*this) - fn(t, 3) * (u.spin() * v.spin() - v.cross(u).spin());
		}
		else
		{
			const double t = norm();
			const int32_t s = q ? -1 : +1;
			const double a = 2 * fn(t, 4) - fn(t, 3);
			const double b = 3 * fn(t, 5) - fn(t, 4);
			return inner(u) * (s * a * spin() + b * spin() * spin()) + s * fn(t, 2) * u.spin() + fn(t, 3) * (u.spin() * spin() + spin() * u.spin());
		}
	}
	Mat3 Vec3::rotation_third_inverse(const Vec3& v, const Vec3& u, bool q, bool x) const
	{
		if(x)
		{
			const double t = norm();
			const double a = (fn(t, 3) - 2 * fn(t, 4)) / fn(t, 2) / 2;
			const double b = (fn(t, 5) - 4 * fn(t, 6)) / fn(t, 2) / 2;
			const double c = (9 * fn(t, 7) - fn(t, 6) - 24 * fn(t, 8)) / fn(t, 2) / 2 + 2 * a * b;
			return
				a * (u.spin() * v.spin() + v.spin() * u.spin()) +
				b * (cross(v).cross(u) - v.cross(cross(u))).outer(*this) - c * cross(v).inner(cross(u)) * outer(*this) +
				b * (outer(cross(u)) * v.spin() + outer(cross(v)) * u.spin() - cross(v).inner(cross(u)) * Mat3::eye());
		}
		else
		{
			const double t = norm();
			const int32_t s = q ? -1 : +1;
			const double a = (fn(t, 3) - 2 * fn(t, 4)) / fn(t, 2) / 2;
			const double b = (fn(t, 5) - 4 * fn(t, 6)) / fn(t, 2) / 2;
			return s * u.spin() / 2 + a * (cross(u).spin() - u.spin() * spin()) - b * outer(cross(u)) * spin();
		}
	}

	//friends
	Vec3 operator*(double s, const Vec3& v)
	{
		return Vec3(v) *= s;
	}
}