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

			//types
			typedef std::function<void(double*, const double*, const double*)> InternalForce;
			typedef std::function<void(double*, const double*, double, double)> ExternalForce;

			typedef std::function<void(double*, const double*)> Inertia;
			typedef std::function<void(double*, const double*, const double*)> Damping;
			typedef std::function<void(double*, const double*, const double*, const double*, double, double, double)> Stiffness;

			//data
			double load(double);
			double load(void) const;

			uint32_t dofs(uint32_t);
			uint32_t dofs(void) const;

			Control control(Control);
			Control control(void) const;

			double frequency(double);
			double frequency(void) const;

			uint32_t harmonics(uint32_t);
			uint32_t harmonics(void) const;

			Inertia inertia(Inertia);
			Damping damping(Damping);
			Stiffness stiffness(Stiffness);

			uint32_t quadrature_order(uint32_t);
			uint32_t quadrature_order(void) const;

			InternalForce internal_force(InternalForce);
			ExternalForce external_force(ExternalForce);

			//tests
			void test_inertia(void) const;
			void test_damping(void) const;
			void test_stiffness(void) const;

			//data
			using Solver::save;
			using Solver::state_set, Solver::force_set, Solver::tangent_set;

			using Solver::silent, Implicit::m_equilibrium;
			using Implicit::convergence, Implicit::continuation, Incremental::stop_criteria;
			using Solver::callback_step, Solver::callback_stop, Solver::callback_record, Solver::callback_update, Solver::callback_restore;

			using Solver::m_watch_dof;
			using Incremental::m_step, Incremental::step_max;
			using Implicit::attempt_max, Implicit::iteration_max;

			using Solver::m_r, Solver::m_fe, Solver::m_K;
			using Solver::m_x_old, Solver::m_x_new, Incremental::m_x_data, Solver::m_dx;
			using Solver::m_p_old, Solver::m_p_new, Incremental::m_p_data, Solver::m_dp, Solver::m_dp0, Solver::m_ddp;

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
			Control m_control;
			Inertia m_inertia;
			Damping m_damping;
			Stiffness m_stiffness;
			InternalForce m_internal_force;
			ExternalForce m_external_force;

			double m_w, m_l;
			double *m_sq, *m_wq;
			double *m_xd, *m_vd, *m_ad;
			double *m_Kd, *m_Cd, *m_Md;
			double *m_rd, *m_fid, *m_fed;
			uint32_t m_dofs, m_harmonics, m_quadrature_order;

		};
	}
}