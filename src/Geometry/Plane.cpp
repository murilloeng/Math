//Math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Linear/Matrix.hpp"
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
			Vector p(4);
			const Vec3& x1 = m_point;
			const Vec3& n1 = m_normal;
			const Vec3& x2 = plane.m_point;
			const Vec3& n2 = plane.m_normal;
			//setup
			Vec3 u1, v1;
			Vec3 u2, v2;
			n1.triad(u1, v1);
			n2.triad(u2, v2);
			//distance
			const Vector b = {
				+(x2 - x1).inner(u1),
				+(x2 - x1).inner(v1),
				-(x2 - x1).inner(u2),
				-(x2 - x1).inner(v2)
			};
			const Matrix A = {
				{1, 0, -u1.inner(u2), -u1.inner(v2)},
				{0, 1, -v1.inner(u2), -v1.inner(v2)},
				{-u2.inner(u1), -u2.inner(v1), 1, 0},
				{-v2.inner(u1), -v2.inner(v1), 0, 1}
			};
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