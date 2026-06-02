#pragma once

//Math
#include "Math/inc/Linear/Vec3.hpp"

namespace math
{
	namespace geometry
	{
		class Circle
		{
		public:
			//constructors
			Circle(void);
			Circle(const Vec3&, const Vec3&, const Vec3&);

			//destructor
			~Circle(void);

			//data
			Vec3 m_center;
			Vec3 m_normal;
			double m_radius;
		};
	}
}