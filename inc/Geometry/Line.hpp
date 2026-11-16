#pragma once

//Math
#include "Math/inc/Linear/Matrix.hpp"

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
			Line(Vec3, Vec3);

			//destructor
			~Line(void) = default;

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