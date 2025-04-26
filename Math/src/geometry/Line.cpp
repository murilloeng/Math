//std
#include <cmath>

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/geometry/Line.hpp"
#include "Math/Math/inc/geometry/Plane.hpp"
#include "Math/Math/inc/geometry/Point.hpp"
#include "Math/Math/inc/geometry/Segment.hpp"

namespace math
{
	namespace geometry
	{
		//constructor
		Line::Line(vec3 point, vec3 direction) : m_point(point), m_direction(direction)
		{
			return;
		}

		//destructor
		Line::~Line(void)
		{
			return;
		}

		//distance
		double Line::distance(const Line& line) const
		{
			//data
			const vec3& x1 = m_point;
			const vec3& t1 = m_direction;
			const vec3& x2 = line.m_point;
			const vec3& t2 = line.m_direction;
			//distance
			const double cm = t1.inner(t2);
			const double r1 = (x2 - x1).inner(t1);
			const double r2 = (x2 - x1).inner(t2);
			const double s1 = +(r1 - cm * r2) / (1 - cm * cm);
			const double s2 = -(r2 - cm * r1) / (1 - cm * cm);
			return (x2 + s2 * t2 - x1 - s1 * t1).norm();
		}
		double Line::distance(const Plane& plane) const
		{
			//data
			const vec3& xl = m_point;
			const vec3& tl = m_direction;
			const vec3& x0 = plane.m_point;
			const vec3& tn = plane.m_normal;
			//distance
			vec3 t1, t2;
			tn.triad(t1, t2);
			const double c1 = tl.inner(t1);
			const double c2 = tl.inner(t2);
			const double r1 = (xl - x0).inner(t1);
			const double r2 = (xl - x0).inner(t2);
			const double rl = (xl - x0).inner(tl);
			const double d = c1 * c1 + c2 * c2 - 1;
			const double s = (rl - c1 * r1 - c2 * r2) / d;
			const double a = ((c2 * c2 - 1) * r1 + c1 * (rl - c2 * r2)) / d;
			const double b = ((c1 * c1 - 1) * r2 + c2 * (rl - c1 * r1)) / d;
			return (x0 + a * t1 + b * t2 - xl - s * tl).norm();
		}
		double Line::distance(const Point& point) const
		{
			return point.distance(*this);
		}
		double Line::distance(const Segment& segment) const
		{
			//data
			const vec3& x1 = m_point;
			const vec3& t1 = m_direction;
			const vec3& x2 = segment.m_points[0];
			const vec3& t2 = segment.direction();
			//distance
			const double cm = t1.inner(t2);
			const double r1 = (x2 - x1).inner(t1);
			const double r2 = (x2 - x1).inner(t2);
			const double s1 = +(r1 - cm * r2) / (1 - cm * cm);
			const double s2 = -(r2 - cm * r1) / (1 - cm * cm);
			return (x2 + bound(s2, 0, 1) * t2 - x1 - s1 * t1).norm();
		}
	}
}