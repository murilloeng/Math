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
				Load, Frequency
			};

			//tests
			void test_inertia(void) const;
			void test_damping(void) const;
			void test_stiffness(void) const;

			//types
			typedef std::function<void(double*, const double*, const double*)> InternalForce;
			typedef std::function<void(double*, const double*, double, double)> ExternalForce;

			typedef std::function<void(double*, const double*)> Inertia;
			typedef std::function<void(double*, const double*, const double*)> Damping;
			typedef std::function<void(double*, const double*, const double*, const double*, double, double, double)> Stiffness;

			//data
			uint32_t dofs(uint32_t);
			uint32_t dofs(void) const;

			uint32_t harmonics(uint32_t);
			uint32_t harmonics(void) const;

			uint32_t quadrature_order(uint32_t);
			uint32_t quadrature_order(void) const;

			double load(double);
			double load(void) const;

			double frequency(double);
			double frequency(void) const;

			Control control(Control);
			Control control(void) const;

			Inertia inertia(Inertia);
			Inertia inertia(void) const;

			Damping damping(Damping);
			Damping damping(void) const;

			Stiffness stiffness(Stiffness);
			Stiffness stiffness(void) const;

			InternalForce internal_force(void) const;
			InternalForce internal_force(InternalForce);

			ExternalForce external_force(void) const;
			ExternalForce external_force(ExternalForce);

			//solve
			void solve(void) override;
			void cleanup(void) override;
			void allocate(void) override;
			void allocate(uint32_t, uint32_t, uint32_t);

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

			//data
			using Solver::size, Solver::allocate;
			using NewtonRaphson::system, NewtonRaphson::residue, NewtonRaphson::tangent_p, NewtonRaphson::tangent_x;

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

			std::function<void(double*, const double*)> m_inertia;
			std::function<void(double*, const double*, const double*)> m_damping;
			std::function<void(double*, const double*, const double*, const double*, double, double, double)> m_stiffness;

		};
	}
}