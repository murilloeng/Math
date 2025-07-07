//std
#include <cmath>

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/groups/ASO3.hpp"
#include "Math/Math/inc/groups/GSO3.hpp"

namespace math
{
	namespace groups
	{
		//constructor
		GSO3::GSO3(void)
		{
			return;
		}
		
		//destructor
		GSO3::~GSO3(void)
		{
			return;
		}

		//logarithm
		ASO3 GSO3::logarithm(void) const
		{
			//angle
			const double t = 2 * acos(bound(m_quaternion[0]));
			//vector
			vec3 aso3;
			const double s = sin(t / 2);
			aso3[0] = s ? t * m_quaternion[1] / s : 0;
			aso3[1] = s ? t * m_quaternion[2] / s : 0;
			aso3[2] = s ? t * m_quaternion[3] / s : 0;
			//return
			return aso3;
		}
	}
}