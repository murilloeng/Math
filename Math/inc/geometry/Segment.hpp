#pragma once

#include "Math/Math/inc/linear/vec3.hpp"

namespace math
{
	namespace geometry
	{
		class Line;
		class Point;

		class Segment
		{
		public:
			//constructors
			Segment(vec3 = {0, 0, 0}, vec3 = {1, 0, 0});

			//destructor
			~Segment(void);

			//length
			double length(void) const;

			//distance
			double distance(const Line&) const;
			double distance(const Point&) const;
			double distance(const Segment&) const;

			//data
			vec3 m_points[2];
		};
	}
}