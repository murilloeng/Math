#pragma once

#include "Math/inc/Linear/Vec3.hpp"

namespace math
{
	namespace geometry
	{
		class Line;
		class Plane;
		class Point;

		class Segment
		{
		public:
			//constructors
			Segment(Vec3 = {0, 0, 0}, Vec3 = {1, 0, 0});

			//destructor
			~Segment(void);

			//length
			double length(void) const;

			//direction
			Vec3 direction(void) const;

			//distance
			double distance(const Line&) const;
			double distance(const Plane&) const;
			double distance(const Point&) const;
			double distance(const Segment&) const;

			//data
			Vec3 m_points[2];
		};
	}
}