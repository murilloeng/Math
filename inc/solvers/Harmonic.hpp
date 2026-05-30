#pragma once

//Math
#include "Math/inc/solvers/NewtonRaphson.hpp"

namespace math
{
	namespace solvers
	{
		class Harmonic : private NewtonRaphson
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

			using NewtonRaphson::state_set;

		private:
			//solve
			void apply(void) override;
			void check(void) override;
			void setup(void) override;
			void compute(void) override;
			void predictor(void) override;
			void corrector(void) override;

			//state
			void compute_state(double);
			void compute_residue(double);
			void compute_velocity(double);
			void compute_acceleration(double);

			//tangent
			void compute_tangent_l(double);
			void compute_tangent_w(double);
			void compute_tangent_z(double);

			//harmonic
			void compute_harmonic_residue(void);
			void compute_harmonic_tangent_p(void);
			void compute_harmonic_tangent_l(void);
			void compute_harmonic_tangent_w(void);
			void compute_harmonic_tangent_z(void);

		public:
			//solver
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
			std::function<void(double*, const double*, const double*, double, double)> m_external_force;

			std::function<void(double*, const double*)> m_inertia;
			std::function<void(double*, const double*, const double*, double)> m_damping;
			std::function<void(double*, const double*, const double*, const double*, double, double, double)> m_stiffness;
		};
	}
}