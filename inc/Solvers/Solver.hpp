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

			//serialization
			virtual void save(const char*) const;

			//data
			virtual uint32_t state_set(void) const = 0;
			virtual uint32_t force_set(void) const = 0;
			virtual uint32_t tangent_set(void) const = 0;

			//enums
			enum class State : uint32_t
			{
				x = 1 << 0, v = 1 << 1, a = 1 << 2, p = 1 << 3, t = 1 << 4
			};
			enum class Force : uint32_t
			{
				r = 1 << 0, fi = 1 << 1, fe = 1 << 2
			};
			enum class Tangent : uint32_t
			{
				K = 1 << 0, C = 1 << 1, M = 1 << 2
			};

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
			//solve
			virtual void step(void);
			virtual void solve(void);
			virtual void cleanup(void);
			virtual void allocate(void);
			virtual void allocate(uint32_t);

			//data
			bool m_silent;
			int32_t* m_rows_map;
			int32_t* m_cols_map;
			uint32_t m_size, m_watch_dof;

			std::function<bool(void)> m_stop;
			std::function<void(void)> m_record;
			std::function<void(void)> m_update;
			std::function<void(void)> m_restore;
			std::function<void(uint32_t)> m_interface;

			double *m_K, *m_C, *m_M;
			double *m_r, *m_fi, *m_fe;
			double *m_dxr, *m_dxt, *m_ddxr, *m_ddxt;

			double *m_x_old, *m_x_new, *m_dx;
			double *m_v_old, *m_v_new, *m_dv;
			double *m_a_old, *m_a_new, *m_da;
			double m_p_old, m_p_new, m_dp, m_dp0, m_ddp;
			double m_t_old, m_t_new, m_dt, m_t_min, m_t_max;
		};
	}
}