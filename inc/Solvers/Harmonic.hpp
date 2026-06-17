#pragma once

//Math
#include "Math/inc/Solvers/NewtonRaphson.hpp"

namespace math
{
	namespace solvers
	{
		class Harmonic : virtual private NewtonRaphson
		{
		public:
			//constructors
			Harmonic(void);

			//destructor
			~Harmonic(void);

			//enums
			enum class Control : uint32_t
			{
				Load,
				Frequency
			};

			//tests
			void test_inertia(void) const;
			void test_damping(void) const;
			void test_stiffness(void) const;

			//data
			using Solver::save, Solver::solve;
			using Solver::state_set, Solver::force_set, Solver::tangent_set;
			
			using Solver::m_silent, Solver::m_equilibrium;
			using Solver::m_convergence, Solver::m_continuation, Solver::m_stop_criteria;
			using Solver::m_stop, Solver::m_record, Solver::m_update, Solver::m_restore, Solver::m_interface;

			using Solver::m_watch_dof;
			using Solver::m_step, Solver::m_attempt, Solver::m_iteration;
			using Solver::m_step_max, Solver::m_attempt_max, Solver::m_iteration_max;

			using Solver::m_r, Solver::m_fe, Solver::m_K;
			using Solver::m_x_old, Solver::m_x_new, Solver::m_x_data, Solver::m_dx;
			using Solver::m_p_old, Solver::m_p_new, Solver::m_p_data, Solver::m_dp, Solver::m_dp0, Solver::m_ddp;

		private:
			//solve
			void apply(void) override;
			void check(void) override;
			void setup(void) override;

			//state
			void compute_state(const double*, double);
			void compute_residue(const double*, double);
			void compute_velocity(const double*, double);
			void compute_acceleration(const double*, double);

			//tangent
			void compute_tangent_l(const double*, double);
			void compute_tangent_w(const double*, double);
			void compute_tangent_z(const double*, double);

			//harmonic
			void compute_harmonic_residue(double*, const double*);
			void compute_harmonic_tangent_p(double*, const double*);
			void compute_harmonic_tangent_l(double*, const double*);
			void compute_harmonic_tangent_w(double*, const double*);
			void compute_harmonic_tangent_z(double*, const double*);

		public:
			//solve
			void solve(void) override;
			void cleanup(void) override;
			void allocate(void) override;

			//data
			double m_w, m_l;
			double *m_sq, *m_wq;
			double *m_xd, *m_vd, *m_ad;
			double *m_Kd, *m_Cd, *m_Md;
			double *m_rd, *m_fid, *m_fed;
			uint32_t m_dofs, m_harmonics, m_quadrature_order;

			Control m_control;
			uint32_t m_stability_steps;
			bool m_stability, *m_stability_data;
			std::function<void(double*, const double*, const double*)> m_internal_force;
			std::function<void(double*, const double*, double, double)> m_external_force;

			std::function<void(double*, const double*)> m_inertia;
			std::function<void(double*, const double*, const double*)> m_damping;
			std::function<void(double*, const double*, const double*, const double*, double, double, double)> m_stiffness;
		};
	}
}