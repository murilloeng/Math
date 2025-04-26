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
		Point::Point(vec3 coordinates) : m_position(coordinates)
		{
			return;
		}
		
		//destructor
		Point::~Point(void)
		{
			return;
		}

		//distance
		double Point::distance(const Line& line) const
		{
			//data
			const vec3& xp = m_position;
			const vec3& xl = line.m_point;
			const vec3& tl = line.m_direction;
			const double a = (xp - xl).inner(tl);
			//distance
			return (xl - xp + a * tl).norm();
		}
		double Point::distance(const Plane& plane) const
		{
			//data
			vec3 t1, t2;
			const vec3& xp = m_position;
			const vec3& x0 = plane.m_point;
			const vec3& tn = plane.m_normal;
			//distance
			tn.triad(t1, t2);
			const double a = (xp - x0).inner(t1);
			const double b = (xp - x0).inner(t2);
			return (x0 - xp + a * t1 + b * t2).norm();
		}
		double Point::distance(const Point& point) const
		{
			return (point.m_position - m_position).norm();
		}
		double Point::distance(const Segment& segment) const
		{
			//data
			const vec3& xp = m_position;
			const vec3& xl = segment.m_points[0];
			const vec3& tl = segment.direction();
			const double a =  (xp - xl).inner(tl);
			//distance
			return (xl - xp + bound(a, 0, 1) * tl).norm();
		}
	}
}