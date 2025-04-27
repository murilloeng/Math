#pragma once

#include "Math/Math/inc/linear/vec3.hpp"

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
			vec3 m_point;
			vec3 m_normal;
		};
	}
}