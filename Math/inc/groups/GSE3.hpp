#pragma once

//math
#include "Math/Math/inc/linear/vec3.hpp"
#include "Math/Math/inc/linear/quat.hpp"

namespace math
{
	namespace groups
	{
		class GSE3
		{
			//constructors
			GSE3(void);
			GSE3(vec3, quat);

			//destructor
			~GSE3(void);

			//matrix
			

		private:
			//data
			vec3 m_vector;
			quat m_quaternion;
		};
	}
}