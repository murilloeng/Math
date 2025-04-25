#pragma once

#include "Math/Math/inc/linear/vec3.hpp"

namespace math
{
	namespace geometry
	{
		class Line;
		class Segment;

		class Point
		{
		public:
			//constructor
			Point(vec3 = {0, 0, 0});

			//destructor
			~Point(void);

			//distance
			double distance(const Line&) const;
			double distance(const Point&) const;
			double distance(const Segment&) const;

			//data
			vec3 m_coordinates;
		};
	}
}