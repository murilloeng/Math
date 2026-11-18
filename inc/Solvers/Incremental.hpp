#pragma once

//Math
#include "Math/inc/Solvers/Solver.hpp"
#include "Math/inc/Solvers/StopCriteria.hpp"

namespace math
{
	namespace solvers
	{
		class Incremental : public virtual Solver
		{
		protected:
			//constructor
			Incremental(void);

			//destructor
			~Incremental(void);

		public:
			//solve
			void solve(void) override;
			void cleanup(void) override;
			virtual void compute_step(void) = 0;

			//data
			uint32_t step(void) const;

			uint32_t step_max(uint32_t);
			uint32_t step_max(void) const;

			bool step_adapt(bool);
			bool step_adapt(void) const;

			double step_scale(void) const;

			double step_scale_min(double);
			double step_scale_min(void) const;

			double step_scale_max(double);
			double step_scale_max(void) const;

			double step_factor_min(double);
			double step_factor_min(void) const;

			double step_factor_max(double);
			double step_factor_max(void) const;

			StopCriteria& stop_criteria(void);

			const double* time_data(void) const;
			const double* load_data(void) const;
			const double* state_data(void) const;
			const double* velocity_data(void) const;
			const double* acceleration_data(void) const;
			const double* kinetic_energy_data(void) const;
			const double* internal_energy_data(void) const;
			const double* potential_energy_data(void) const;

			double time_data(uint32_t) const;
			double load_data(uint32_t) const;

			double kinetic_energy_data(uint32_t) const;
			double internal_energy_data(uint32_t) const;
			double potential_energy_data(uint32_t) const;

			double state_data(uint32_t, uint32_t) const;
			double velocity_data(uint32_t, uint32_t) const;
			double acceleration_data(uint32_t, uint32_t) const;

			bool stability(bool);
			bool stability(void) const;

			bool stability_data(uint32_t) const;
			const bool* stability_data(void) const;

			//serialization
			void save(const char*) const override;
			void save(const char*, std::vector<uint32_t>) const;

		protected:
			//solve
			bool stop(void) override;
			void print(void) override;
			void setup(void) override;
			void record(void) override;

			//allocate
			void allocate_state(void) override;
			void allocate_energies(void) override;

			//data
			bool m_step_adapt, m_step_stop;
			uint32_t m_step, m_step_max;
			StopCriteria m_stop_criteria;
			bool m_stability, *m_stability_data;
			double m_step_factor_min, m_step_factor_max;
			double m_step_scale, m_step_scale_min, m_step_scale_max;

			double *m_p_data, *m_t_data;
			double *m_x_data, *m_v_data, *m_a_data;
			double *m_T_data, *m_U_data, *m_V_data;

			//friends
			friend class math::solvers::StopCriteria;
		};
	}
}