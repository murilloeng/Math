#pragma once

namespace tests
{
	namespace fem
	{
		namespace beam
		{
			namespace dynamics
			{
				void section_strains_gradient(void);
			}
		}
		void revolute_fixed(void);
		void revolute_flexible(void);
	}
	namespace misc
	{
		void fft(void);
		void drift(void);
	}
	namespace groups
	{
		namespace GSO3
		{
			void log(void);
			void inverse(void);
			void tangent(void);
			void tangent_inverse(void);
			void tangent_indentity(void);
			void tangent_increment(void);
			void tangent_inverse_increment(void);
		}
		namespace GSE3
		{
			void log(void);
			void inverse(void);
			void tangent(void);
			void tangent_inverse(void);
			void tangent_increment(void);
			void tangent_inverse_increment(void);
		}
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
		namespace harmonic
		{
			void pyramid(void);
			void oscillator(void);
		}
		namespace newton_raphson
		{
			void truss_von_mises(void);
		}
		namespace newmark
		{
			void single_dof(void);
			void single_pendulum(void);
		}
	}
	namespace geometry
	{
		void circle(void);
	}
}