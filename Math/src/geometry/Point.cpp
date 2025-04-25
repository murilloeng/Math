//math
#include "Math/Math/inc/geometry/Point.hpp"

namespace math
{
	namespace geometry
	{
		//constructor
		Point::Point(vec3 coordinates) : m_coordinates(coordinates)
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
			return 0;
		}
		double Point::distance(const Point& point) const
		{
			return (point.m_coordinates - m_coordinates).norm();
		}
		double Point::distance(const Segment& segment) const
		{
			return 0;
		}
	}
}