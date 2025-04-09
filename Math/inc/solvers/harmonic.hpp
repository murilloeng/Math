#pragma once

//math
#include "Math/Math/inc/quadrature/rule.hpp"
#include "Math/Math/inc/quadrature/quadrature.hpp"

namespace math
{
	class harmonic
	{
	public:
		//constructors
		harmonic(void);

		//destructor
		~harmonic(void);

		//solve
		bool solve(void);
		void setup(void);
		void cleanup(void);

		//data
		const double* amplitudes(void) const;

		//data
		void** m_args;
		uint32_t m_size;
		uint32_t m_harmonics;
		uint32_t m_iteration_max;
		uint32_t m_quadrature_order;
		double m_frequency, m_tolerance;

		void(*m_external_force)(double*, double, const double*, void**);
		void(*m_internal_force)(double*, const double*, const double*, void**);

		void(*m_inertia)(double*, const double*, void**);
		void(*m_damping)(double*, const double*, const double*, void**);
		void(*m_stiffness)(double*, double, const double*, const double*, void**);

	private:
		//solver
		void compute_residue(void);
		void compute_tangent(void);

		void compute_state(double);
		void compute_residue(double);
		void compute_velocity(double);
		void compute_acceleration(double);

		//data
		quadrature::Quadrature* m_quadrature;
		double *m_d, *m_v, *m_a, *m_z, *m_r, *m_fi, *m_fe, *m_fr, *m_Kt, *m_Ct, *m_Mt, *m_At;
	};
}