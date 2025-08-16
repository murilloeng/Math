#pragma once

namespace tests
{
	namespace rotations
	{
		namespace vec3
		{
			void rotation_tensor(void);
			void rotation_hessian(void);
			void rotation_gradient(void);
		}
		namespace quat
		{
			void rotation_tensor(void);
		}
	}
}