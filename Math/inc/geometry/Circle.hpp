#pragma once

//Math
#include "Math/Math/inc/linear/vec3.hpp"

namespace math
{
	namespace geometry
	{
		class Circle
		{
		public:
			//constructors
			Circle(void);
			Circle(const vec3&, const vec3&, const vec3&);

			//destructor
			~Circle(void);

			//data
			vec3 m_center;
			vec3 m_normal;
			double m_radius;
		};
	}
}