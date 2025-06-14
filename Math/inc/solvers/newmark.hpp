#pragma once

//std
#include <cstdint>

//Math
#include "Math/Math/inc/solvers/solver.hpp"

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
		class newmark : public solver
		{
		public:
			//constructors
			newmark(void);

			//destructor
			~newmark(void);

			//data
			uint32_t state_set(void) const override;
			uint32_t force_set(void) const override;
			uint32_t tangent_set(void) const override;

		private:
			//solve
			void check(void) override;
			void setup(void) override;
			void compute(void) override;
			void predictor(void) override;
			void corrector(void) override;

		public:
			//data
			double m_g, m_b;

			std::function<void(double*, const double*, const double*)> m_internal_force;
			std::function<void(double*, const double*, const double*, double)> m_external_force;

			std::function<void(double*, const double*)> m_inertia;
			std::function<void(double*, const double*, const double*, double)> m_damping;
			std::function<void(double*, const double*, const double*, const double*, double)> m_stiffness;
		};
	}
}