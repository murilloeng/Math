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
				x = 1 << 0,
				v = 1 << 1,
				a = 1 << 2,
				p = 1 << 3,
				t = 1 << 4
			};
			enum class Force : uint32_t
			{
				r = 1 << 0,
				fi = 1 << 1,
				fe = 1 << 2
			};
			enum class Tangent : uint32_t
			{
				K = 1 << 0,
				C = 1 << 1,
				M = 1 << 2
			};

			//data
			bool silent(bool);
			bool silent(void) const;

			uint32_t size(uint32_t);
			uint32_t size(void) const;
			
			int32_t* rows_map(int32_t*);
			int32_t* rows_map(void) const;
			
			int32_t* cols_map(int32_t*);
			int32_t* cols_map(void) const;
			
			uint32_t watch_dof(uint32_t);
			uint32_t watch_dof(void) const;

			double time_min(double);
			double time_min(void) const;

			double time_max(double);
			double time_max(void) const;

			double load_increment(double);
			double load_increment(void) const;
			
			double* state_old(void) const;
			double* state_old(const double*);
			
			double* state_new(void) const;
			double* state_new(const double*);

			double* velocity_old(void) const;
			double* velocity_old(const double*);

			double* velocity_new(void) const;
			double* velocity_new(const double*);

			double* acceleration_old(void) const;
			double* acceleration_old(const double*);

			double* acceleration_new(void) const;
			double* acceleration_new(const double*);

			std::function<void(void)> callback_step(std::function<void(void)>);
			std::function<bool(void)> callback_stop(std::function<bool(void)>);
			std::function<void(void)> callback_record(std::function<void(void)>);
			std::function<void(void)> callback_update(std::function<void(void)>);
			std::function<void(void)> callback_restore(std::function<void(void)>);

			//serialization
			virtual void save(const char*) const;

			//sets
			virtual uint32_t state_set(void) const = 0;
			virtual uint32_t force_set(void) const = 0;
			virtual uint32_t tangent_set(void) const = 0;

			//solve
			virtual void step(void);
			virtual void solve(void);
			virtual void cleanup(void);
			virtual void allocate(void);
			virtual void allocate(uint32_t);

		protected:
			//solve
			virtual bool stop(void);
			virtual void apply(void);
			virtual void print(void);
			virtual void setup(void);
			virtual void record(void);
			virtual void update(void);
			virtual void restore(void);

			//solve
			virtual void check(void) = 0;
			virtual void compute(void) = 0;
			virtual void predictor(void) = 0;
			virtual void corrector(void) = 0;

			//allocate
			virtual void allocate_state(void);
			virtual void allocate_forces(void);
			virtual void allocate_tangents(void);

			//solve
			bool solve(const double*, const double*, double*) const;

			//data
			bool m_silent;
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
			double *m_dxr, *m_dxt, *m_ddxr, *m_ddxt;

			double *m_x_old, *m_x_new, *m_dx;
			double *m_v_old, *m_v_new, *m_dv;
			double *m_a_old, *m_a_new, *m_da;
			double m_p_old, m_p_new, m_dp, m_dp0, m_ddp;
			double m_t_old, m_t_new, m_dt, m_t_min, m_t_max;

			//friends
			friend class Convergence;
			friend class Continuation;
		};
	}
}