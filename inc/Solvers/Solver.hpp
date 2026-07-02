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

			//sets
			virtual uint32_t state_set(void) const = 0;
			virtual uint32_t force_set(void) const = 0;
			virtual uint32_t tangent_set(void) const = 0;

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

		public:
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
			double *m_x_old, *m_x_new, *m_dx;
			double *m_v_old, *m_v_new, *m_dv;
			double *m_a_old, *m_a_new, *m_da;
			double *m_dxr, *m_dxt, *m_ddxr, *m_ddxt;
			double m_p_old, m_p_new, m_dp, m_dp0, m_ddp;
			double m_t_old, m_t_new, m_dt, m_t_min, m_t_max;
		};
	}
}