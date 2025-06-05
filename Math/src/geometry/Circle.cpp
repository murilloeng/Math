//std
#include <cmath>

//Math
#include "Math/Math/inc/geometry/Circle.hpp"

namespace math
{
	namespace geometry
	{
		//constructor
		Circle::Circle(void) : m_center(0, 0, 0), m_normal(0, 0, 1), m_radius(1)
		{
			return;
		}
		Circle::Circle(const vec3& x1, const vec3& x2, const vec3& x3)
		{
			//normal
			m_normal = (x2 - x1).cross(x3 - x1).normalize();
			//center
			matrix A(3, 3), f(3, 1);
			for(uint32_t i = 0; i < 3; i++)
			{
				A(2, i) = m_normal[i];
				A(0, i) = 2 * (x2[i] - x1[i]);
				A(1, i) = 2 * (x3[i] - x1[i]);
				f(2) = m_normal.inner(x1);
				f(0) = x2.inner(x2) - x1.inner(x1);
				f(1) = x3.inner(x3) - x1.inner(x1);
			}
			A.solve(m_center, f);
			//radius
			m_radius = (x1 - m_center).norm();
		}
		
		//destructor
		Circle::~Circle(void)
		{
			return;
		}
	}
}