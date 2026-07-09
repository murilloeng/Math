#pragma once

//std
#include <cstdint>
#include <functional>

//Math
#include "Math/inc/Solvers/StopCriteria.hpp"

namespace math
{
	namespace solvers
	{
		class Convergence;
		class Continuation;
	}
}

namespace math
{
	namespace solvers
	{
		class Solver
		{
		public:
			//constructor
			Solver(void);

			//destructor
			virtual ~Solver(void);

			//enums
			enum class State : uint32_t
			{
				x, v, a, p, t
			};
			enum class Force : uint32_t
			{
				r, fi, fe
			};
			enum class Tangent : uint32_t
			{
				K, C, M
			};

			//solve
			virtual void solve(void) = 0;

			//data
			virtual void cleanup(void);
			virtual void allocate(void);
			virtual void allocate(uint32_t);

			//serialization
			virtual void save(const char*) const;

			//data
			virtual uint32_t state_set(void) const = 0;
			virtual uint32_t force_set(void) const = 0;
			virtual uint32_t tangent_set(void) const = 0;

			//data
			bool silent(bool);
			bool silent(void) const;

			bool status(void) const;

			uint32_t size(uint32_t);
			uint32_t size(void) const;

			uint32_t watch_dof(uint32_t);
			uint32_t watch_dof(void) const;

			int32_t* rows_map(int32_t*);
			int32_t* cols_map(int32_t*);
			int32_t* rows_map(void) const;
			int32_t* cols_map(void) const;

			double time_min(double);
			double time_max(double);
			double time_min(void) const;
			double time_max(void) const;

			double step_size(double);
			double step_size(void) const;

			double* state_old(void);
			double* state_new(void);
			double state_old(uint32_t);
			double state_new(uint32_t);
			double state_old(uint32_t, double);
			double state_new(uint32_t, double);

			double* velocity_old(void);
			double* velocity_new(void);
			double velocity_old(uint32_t);
			double velocity_new(uint32_t);
			double velocity_old(uint32_t, double);
			double velocity_new(uint32_t, double);

			double* acceleration_old(void);
			double* acceleration_new(void);
			double acceleration_old(uint32_t);
			double acceleration_new(uint32_t);
			double acceleration_old(uint32_t, double);
			double acceleration_new(uint32_t, double);

		protected:
			//solve
			virtual bool stop(void);
			virtual void check(void);
			virtual void apply(void);
			virtual void print(void);
			virtual void setup(void);
			virtual void record(void);
			virtual void update(void);
			virtual void restore(void);
			virtual void compute(void);
			virtual void predictor(void);
			virtual void corrector(void);

			//allocate
			virtual void allocate_state(void);
			virtual void allocate_forces(void);
			virtual void allocate_tangents(void);

			//solve
			bool solve(const double*, const double*, double*) const;

			//data
			bool m_silent;
			bool m_status;
			int32_t* m_rows_map;
			int32_t* m_cols_map;
			uint32_t m_size, m_watch_dof;

			std::function<void(void)> m_callback_step;
			std::function<bool(void)> m_callback_stop;
			std::function<void(void)> m_callback_record;
			std::function<void(void)> m_callback_update;
			std::function<void(void)> m_callback_restore;

			double *m_K, *m_C, *m_M;
			double *m_r, *m_fi, *m_fe;
			double *m_x_old, *m_x_new, *m_dx;
			double *m_v_old, *m_v_new, *m_dv;
			double *m_a_old, *m_a_new, *m_da;
			double *m_dxr, *m_dxt, *m_ddxr, *m_ddxt;
			double m_p_old, m_p_new, m_dp, m_dp0, m_ddp;
			double m_t_old, m_t_new, m_dt, m_t_min, m_t_max;

			//friends
			friend class math::solvers::Convergence;
			friend class math::solvers::Continuation;
		};
	}
}