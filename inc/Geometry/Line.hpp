#pragma once

#include "Math/inc/Linear/Vec3.hpp"

namespace math
{
	namespace geometry
	{
		class Point;
		class Plane;
		class Segment;

		class Line
		{
		public:
			//constructors
			Line(Vec3 = {0, 0, 0}, Vec3 = {1, 0, 0});

			//destructor
			~Line(void);

			//distance
			double distance(const Line&) const;
			double distance(const Plane&) const;
			double distance(const Point&) const;
			double distance(const Segment&) const;

			//data
			Vec3 m_point;
			Vec3 m_direction;
		};
	}
}