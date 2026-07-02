#pragma once

//Math
#include "Math/inc/Solvers/NewtonRaphson.hpp"

namespace math
{
	namespace solvers
	{
		class Harmonic : virtual public NewtonRaphson
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
			typedef std::function<void(double*, const double*, const double*, const double*, double, double, double)> Stiffness;

			//tests
			void test_inertia(void) const;
			void test_damping(void) const;
			void test_stiffness(void) const;

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

			//derived
			using Solver::m_size;
			using Solver::allocate;
			using NewtonRaphson::m_residue;
			using NewtonRaphson::m_system_1, NewtonRaphson::m_system_2;
			using NewtonRaphson::m_tangent_1, NewtonRaphson::m_tangent_2;

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
			std::function<void(double*, const double*, const double*)> m_internal_force;
			std::function<void(double*, const double*, double, double)> m_external_force;

			std::function<void(double*, const double*)>  m_inertia;
			std::function<void(double*, const double*, const double*)> m_damping;
			std::function<void(double*, const double*, const double*, const double*, double, double, double)> m_stiffness;

		};
	}
}