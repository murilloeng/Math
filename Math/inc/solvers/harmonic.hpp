#pragma once

//math
#include "Math/Math/inc/quadrature/rule.hpp"

namespace math
{
	enum class harmonic_control
	{
		load		= 0,
		frequency	= 1
	};
	enum class harmonic_strategy
	{
		minimal_norm			= 0,
		uniform_increment		= 1,
		arc_length_spherical	= 2,
		arc_length_cylindrical	= 3
	};
}

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

		//save
		void save(void) const;

		//test
		void test_damping(uint32_t) const;
		void test_stiffness(uint32_t) const;

		//data
		void** m_args;
		uint32_t m_size;
		uint32_t m_step_max;
		uint32_t m_harmonics;
		uint32_t m_attempt_max;
		uint32_t m_iteration_max;
		uint32_t m_amplitude_steps;
		uint32_t m_stability_steps;
		uint32_t m_quadrature_order;

		bool m_stability;
		harmonic_control m_control;
		harmonic_strategy m_strategy;

		double m_dpg, m_tolerance;
		double m_l_0, m_l_min, m_l_max;
		double m_w_0, m_w_min, m_w_max;

		void(*m_internal_force)(double*, const double*, const double*, void**);
		void(*m_external_force)(double*, double, double, const double*, void**);

		void(*m_inertia)(double*, const double*, void**);
		void(*m_damping)(double*, const double*, const double*, void**);
		void(*m_stiffness)(double*, double, double, double, const double*, const double*, const double*, void**);

	private:
		//solver
		bool stop(void);
		void setup(void);
		void cleanup(void);
		void allocate(void);
		void initialize(void);

		//state
		void compute_state(double);
		void compute_residue(void);
		void compute_velocity(double);
		void compute_acceleration(double);

		//system
		void compute_stability(void);
		void compute_system_residue(void);
		void compute_system_tangent_p(void);
		void compute_system_tangent_l(void);
		void compute_system_tangent_w(void);
		void compute_system_tangent_z(void);

		//test
		static void function_residue_v(double*, const double*, const void**);
		static void function_residue_d(double*, const double*, const void**);

		//solver
		void update(void);
		void restore(void);
		void record(uint32_t);
		
		void increment_state(void);
		void compute_predictor(void);
		void compute_corrector(void);
		void compute_amplitudes(void);

		void compute_parameter_predictor(void);
		void compute_parameter_corrector(void);

		void predictor_minimal_norm(void);
		void corrector_minimal_norm(void);

		void predictor_uniform_increment(void);
		void corrector_uniform_increment(void);

		void predictor_arc_length_spherical(void);
		void corrector_arc_length_spherical(void);

		void predictor_arc_length_cylindrical(void);
		void corrector_arc_length_cylindrical(void);

		//data
		uint32_t m_step;
		uint32_t m_attempt;
		uint32_t m_iteration;
		bool m_system_equilibrium;

		double *m_sq, *m_wq;
		double *m_d, *m_v, *m_a;
		
		double m_dp, m_ddp;
		double m_l_old, m_l_new, *m_l_data;
		double m_w_old, m_w_new, *m_w_data;
		double *m_z_old, *m_z_new, *m_z_data;
		double *m_dz, *m_dz0r, *m_dz0t, *m_ddzr, *m_ddzt;

		bool* m_stability_data;
		double *m_amplitudes_data;
		double *m_Kt, *m_Ct, *m_Mt, *m_At, *m_bt;
		double *m_r, *m_fi, *m_fe, *m_fr, *m_dfrw;
		double *m_sd, *m_sv, *m_sa, *m_sr, *m_sc, *m_St, *m_Sm, *m_Sz;
	};
}