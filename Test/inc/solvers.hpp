#pragma once

namespace tests
{
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
		namespace gradient_descent
		{
			void single_quartic(void);
			void single_quadratic(void);
			void double_quadratic(void);
			void exponential_smooth(void);
			void rosenbrock_function(void);
			void himmelblau_function(void);
		}
	}
}