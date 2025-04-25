//math
#include "Math/Math/inc/geometry/Line.hpp"

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
			return 0;
		}
		double Line::distance(const Point& point) const
		{
			return 0;
		}
		double Line::distance(const Segment& segment) const
		{
			return 0;
		}
	}
}