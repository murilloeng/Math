#pragma once

namespace tests
{
	namespace fem
	{
		void revolute_fixed(void);
		void revolute_flexible(void);
	}
	namespace misc
	{
		void fft(void);
		void drift(void);
	}
	namespace rotations
	{
		void vec3_rotation_tensor(void);
		void quat_rotation_tensor(void);
		void vec3_rotation_hessian(void);
		void vec3_rotation_gradient(void);
	}
	namespace solvers
	{
		void harmonic_pyramid(void);
		void harmonic_oscillator(void);
	}
}