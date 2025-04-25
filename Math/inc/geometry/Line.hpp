#pragma once

#include "Math/Math/inc/linear/vec3.hpp"

namespace math
{
	namespace geometry
	{
		class Point;
		class Segment;

		class Line
		{
		public:
			//constructors
			Line(vec3 = {0, 0, 0}, vec3 = {1, 0, 0});

			//destructor
			~Line(void);

			//distance
			double distance(const Line&) const;
			double distance(const Point&) const;
			double distance(const Segment&) const;

			//data
			vec3 m_point;
			vec3 m_direction;
		};
	}
}