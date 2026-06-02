#pragma once

#include "Math/inc/Linear/Vec3.hpp"

namespace math
{
	namespace geometry
	{
		class Line;
		class Point;
		class Segment;

		class Plane
		{
		public:
			//constructor
			Plane(void);

			//destructor
			~Plane(void);

			//distance
			double distance(const Line&) const;
			double distance(const Plane&) const;
			double distance(const Point&) const;
			double distance(const Segment&) const;

			//data
			Vec3 m_point;
			Vec3 m_normal;
		};
	}
}