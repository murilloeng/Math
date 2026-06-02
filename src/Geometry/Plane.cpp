//std
#include <cmath>

//Math
#include "Math/inc/Geometry/Line.hpp"
#include "Math/inc/Geometry/Point.hpp"
#include "Math/inc/Geometry/Plane.hpp"
#include "Math/inc/Geometry/Segment.hpp"

namespace math
{
	namespace geometry
	{
		//constructor
		Plane::Plane(void) : m_point{0, 0, 0}, m_normal{0, 0, 1}
		{
			return;
		}

		//destructor
		Plane::~Plane(void)
		{
			return;
		}

		//distance
		double Plane::distance(const Line& line) const
		{
			return line.distance(*this);
		}
		double Plane::distance(const Point& point) const
		{
			return point.distance(*this);
		}
		double Plane::distance(const Plane& plane) const
		{
			//data
			const Vec3& x1 = m_point;
			const Vec3& n1 = m_normal;
			const Vec3& x2 = plane.m_point;
			const Vec3& n2 = plane.m_normal;
			//setup
			Vec3 u1, v1;
			Vec3 u2, v2;
			n1.triad(u1, v1);
			n2.triad(u2, v2);
			Vector b(4), p(4);
			Matrix A(4, 4, mode::eye);
			//distance
			b(0) = +(x2 - x1).inner(u1);
			b(1) = +(x2 - x1).inner(v1);
			b(2) = -(x2 - x1).inner(u2);
			b(3) = -(x2 - x1).inner(v2);
			A(0, 2) = A(2, 0) = -u1.inner(u2);
			A(0, 3) = A(3, 0) = -u1.inner(v2);
			A(1, 3) = A(3, 1) = -v1.inner(u2);
			A(1, 4) = A(4, 1) = -v1.inner(v2);
			//distance
			A.solve(p, b);
			return (x2 + p[2] * u2 + p[3] * v2 - x1 - p[0] * u1 - p[1] * v1).norm();
		}
		double Plane::distance(const Segment& segment) const
		{
			return segment.distance(*this);
		}
	}
}