#pragma once

//std
#include <cstdint>

//Math
#include "Math/Math/inc/solvers/convergence.hpp"
#include "Math/Math/inc/solvers/stop_criteria.hpp"

//x: state vector
//r: residue vector
//v: velocity vector
//a: acceleration vector
//p: continuation parameter

//target system: r(x(t), v(t), a(t)) = fe(x(t), v(t), t) - fi(x(t), v(t)) - M(x(t)) * a(t) = 0

//tangent on a: M(x) = -dr/da(x)
//tangent on v: C(x, v, t) = -dr/dv(x, v, t) = dfi/dv(x, v) - dfe/dv(x, v, t)
//tangent on x: K(x, v, a) = -dr/dx(x, v, a, t) = dfi/dx(x, v) - dfe/dx(x, v, t) + dM/dx(x) : a

//time interpolation:
//v_new = v_old + dt * a_old + g * dt * (a_new - a_old)
//x_new = x_old + dt * v_old + dt * dt / 2 * a_old + b * dt * dt * (a_new - a_old)

namespace math
{
	namespace solvers
	{
		class newmark
		{
		public:
			//constructors
			newmark(void);

			//destructor
			~newmark(void);

		private:
			//data
			// void update(void);
			void residue(void);
			void inertia(void);
			void damping(void);
			void stifness(void);
			void internal(void);
			void external(void);

			//solve
			bool stop(void);
			void check(void);
			void apply(void);
			void print(void);
			void setup(void);
			void record(void);
			void update(void);
			void restore(void);
			void compute(void);
			void predictor(void);
			void corrector(void);
			bool equilibrium(void);

		public:
			//solve
			void step(void);
			void solve(void);
			void cleanup(void);
			void allocate(void);

			//data
			void** m_args;
			bool m_silent;
			convergence m_convergence;
			stop_criteria m_stop_criteria;

			uint32_t m_watch_dof, m_size;
			uint32_t m_step, m_attempt, m_iteration;
			uint32_t m_step_max, m_attempt_max, m_iteration_max;

			double m_t, m_T, m_dt, m_g, m_b;
			double *m_x_old, *m_x_new, *m_x_data, *m_dx;
			double *m_v_old, *m_v_new, *m_v_data, *m_dv;
			double *m_a_old, *m_a_new, *m_a_data, *m_da;
			double *m_r, *m_fe, *m_fi, *m_M, *m_C, *m_K;

			void (*m_internal_force)(double*, const double*, const double*, void**);
			void (*m_external_force)(double*, const double*, const double*, double, void**);

			void (*m_inertia)(double*, const double*, void**);
			void (*m_damping)(double*, const double*, const double*, double, void**);
			void (*m_stiffness)(double*, const double*, const double*, const double*, double, void**);
		};
	}
}