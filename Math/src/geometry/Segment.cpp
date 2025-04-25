//math
#include "Math/Math/inc/geometry/Line.hpp"
#include "Math/Math/inc/geometry/Point.hpp"
#include "Math/Math/inc/geometry/Segment.hpp"

namespace math
{
	namespace geometry
	{
		//constructor
		Segment::Segment(vec3 point_1, vec3 point_2) : m_points{point_1, point_2}
		{
			return;
		}

		//destructor
		Segment::~Segment(void)
		{
			return;
		}

		//length
		double Segment::length(void) const
		{
			return (m_points[1] - m_points[0]).norm();
		}

		//distance
		double Segment::distance(const Line& line) const
		{
			return line.distance(*this);
		}
		double Segment::distance(const Point& point) const
		{
			return point.distance(*this);
		}
		double Segment::distance(const Segment& segment) const
		{
			return 0;
		}
	}
}