#pragma once

#include "Math/inc/Linear/Vec3.hpp"

namespace math
{
	namespace geometry
	{
		class Line;
		class Plane;
		class Segment;

		class Point
		{
		public:
			//constructor
			Point(Vec3 = {0, 0, 0});

			//destructor
			~Point(void);

			//distance
			double distance(const Line&) const;
			double distance(const Plane&) const;
			double distance(const Point&) const;
			double distance(const Segment&) const;

			//data
			Vec3 m_position;
		};
	}
}