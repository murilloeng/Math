// //std
// #include <cmath>
// #include <cstdio>

// //math
// #include "Math/Math/inc/solvers/harmonic.hpp"

// //test
// #include "Math/Test/inc/tests.hpp"

// static const double k = 4.00e+00;
// static const double c = 1.00e-01;
// static const double m = 1.00e+00;

// static void inertia(double* M, const double* d, void** args)
// {
// 	M[0] = m;
// }
// static void damping(double* C, const double* d, const double* v, void** args)
// {
// 	C[0] = c;
// }
// static void stiffness(double* K, double t, const double* d, const double* v, void** args)
// {
// 	const double z = 1 + d[0];
// 	K[0] = k / 2 * (3 * z * z - 1);
// }

// static void external_force(double* fe, double t, double, const double* d, void** args)
// {
// 	//data
// 	const double w = ((math::harmonic*) args[0])->m_frequency;
// 	//force
// 	fe[0] = 0.04 * k * sqrt(3) / 9 * cos(w * t);
// }
// static void internal_force(double* fi, const double* d, const double* v, void** args)
// {
// 	const double z = 1 + d[0];
// 	fi[0] = c * v[0] + k / 2 * z * (z * z - 1);
// }

// void tests::solvers::harmonic_pyramid(void)
// {
// 	//data
// 	math::harmonic solver;
// 	const uint32_t nh = 5;
// 	const uint32_t nw = 1000;
// 	void* args[] = { &solver };
// 	FILE* file = fopen("harmonic.txt", "w");
// 	//setup
// 	solver.m_size = 1;
// 	solver.m_args = args;
// 	solver.m_harmonics = nh;
// 	solver.m_tolerance = 1e-3;
// 	solver.m_iteration_max = 40;
// 	solver.m_quadrature_order = 200;
// 	//system
// 	solver.m_inertia = inertia;
// 	solver.m_damping = damping;
// 	solver.m_stiffness = stiffness;
// 	solver.m_internal_force = internal_force;
// 	solver.m_external_force = external_force;
// 	//solve
// 	for(uint32_t i = 0; i < nw; i++)
// 	{
// 		//setup
// 		solver.m_frequency = (0.1 + 1.8 * i / (nw - 1)) * sqrt(k / m);
// 		//solve
// 		if(!solver.solve())
// 		{
// 			printf("Harmonic oscillator failed to solve\n");
// 			break;
// 		}
// 		//write
// 		const double* z = solver.amplitudes();
// 		for(uint32_t j = 0; j < 2 * nh + 1; j++)
// 		{
// 			fprintf(file, "%+.6e ", z[j]);
// 		}
// 		fprintf(file, "\n");
// 	}
// 	//close
// 	fclose(file);
// }