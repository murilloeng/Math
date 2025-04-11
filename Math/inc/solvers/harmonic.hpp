#pragma once

//math
#include "Math/Math/inc/quadrature/rule.hpp"

namespace math
{
	enum class harmonic_parameter
	{
		load = 0,
		frequency = 1
	};
	enum class harmonic_strategy
	{
		uniform_increment = 0,
		arc_length_spheric = 1,
		arc_length_cylindric = 2
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
		void save(void);
		bool solve(void);

		//data
		void** m_args;
		uint32_t m_size;
		uint32_t m_step_max;
		uint32_t m_harmonics;
		uint32_t m_attempt_max;
		uint32_t m_iteration_max;
		uint32_t m_quadrature_order;
		harmonic_strategy m_strategy;
		harmonic_parameter m_parameter;
		double m_dpg, m_l_0, m_w_0, m_tolerance;

		void(*m_internal_force)(double*, const double*, const double*, void**);
		void(*m_external_force)(double*, double, double, const double*, void**);
		void(*m_external_force_gradient)(double*, double, double, const double*, void**);

		void(*m_inertia)(double*, const double*, void**);
		void(*m_damping)(double*, const double*, const double*, void**);
		void(*m_stiffness)(double*, double, double, double, const double*, const double*, const double*, void**);

	private:
		//solver
		void setup(void);
		void cleanup(void);
		void allocate(void);
		void initialize(void);

		//state
		void compute_jerk(double);
		void compute_state(double);
		void compute_residue(void);
		void compute_velocity(double);
		void compute_acceleration(double);
		void compute_frequency_gradients(double);

		//system
		void compute_system_residue(void);
		void compute_system_tangent_p(void);
		void compute_system_tangent_l(void);
		void compute_system_tangent_w(void);
		void compute_system_tangent_z(void);

		//solver
		void update(void);
		void restore(void);
		void record(uint32_t);

		void increment_state(void);
		void compute_predictor(void);
		void compute_corrector(void);

		void compute_parameter_predictor(void);
		void compute_parameter_corrector(void);

		void predictor_uniform_increment(void);
		void corrector_uniform_increment(void);

		void predictor_arc_length_spheric(void);
		void corrector_arc_length_spheric(void);

		void predictor_arc_length_cylindric(void);
		void corrector_arc_length_cylindric(void);

		//data
		uint32_t m_step;
		uint32_t m_attempt;
		uint32_t m_iteration;
		bool m_system_equilibrium;

		double *m_sq, *m_wq;
		double *m_d, *m_v, *m_a, *m_j;
		double *m_ddw, *m_dvw, *m_daw;

		double m_l_old, m_l_new, *m_l_data;
		double m_w_old, m_w_new, *m_w_data;
		double *m_z_old, *m_z_new, *m_z_data;

		double m_dp, m_ddp;
		double *m_dz, *m_dz0r, *m_dz0t, *m_ddzr, *m_ddzt;

		double *m_r, *m_fi, *m_fe, *m_fr;
		double *m_Kt, *m_Ct, *m_Mt, *m_At, *m_bt, *m_dfew, *m_dfrw;
	};
}