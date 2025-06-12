#pragma once

//std
#include <functional>

//math
#include "Math/Math/inc/solvers/convergence.hpp"
#include "Math/Math/inc/solvers/continuation.hpp"
#include "Math/Math/inc/solvers/stop_criteria.hpp"

namespace math
{
	namespace solvers
	{
		class nonlinear
		{
		public:
			//constructor
			nonlinear(void);

			//destructor
			virtual ~nonlinear(void);

			//data
			virtual void save(const char*) const = 0;
			virtual uint32_t state_set(void) const = 0;
			virtual uint32_t force_set(void) const = 0;
			virtual uint32_t tangent_set(void) const = 0;

			//enums
			enum class state : uint32_t
			{
				x = 1 << 0, v = 1 << 1, a = 1 << 2, p = 1 << 3
			};
			enum class force : uint32_t
			{
				r = 1 << 0, fi = 1 << 1, fe = 1 << 2
			};
			enum class tangent : uint32_t
			{
				K = 1 << 0, C = 1 << 1, M = 1 << 2
			};

		private:
			//solve
			bool stop(void);
			void check(void);
			void apply(void);
			void print(void);
			void setup(void);
			void record(void);
			void update(void);
			void restore(void);
			void compute(void);
			void predictor(void);
			void corrector(void);
			bool equilibrium(void);

		public:
			//solve
			virtual void step(void);
			virtual void solve(void);
			virtual void cleanup(void);
			virtual void allocate(void);

			//data
			bool m_silent;
			bool m_equilibrium;

			convergence m_convergence;
			continuation m_continuation;
			stop_criteria m_stop_criteria;

			std::function<bool(void)> m_stop;
			std::function<void(void)> m_record;
			std::function<void(void)> m_update;
			std::function<void(void)> m_restore;
			std::function<void(uint32_t)> m_interface;

			uint32_t m_size, m_watch_dof;
			uint32_t m_step, m_attempt, m_iteration;
			uint32_t m_step_max, m_attempt_max, m_iteration_max;

			double *m_K, *m_C, *m_M;
			double *m_r, *m_fi, *m_fe;
			double *m_dxr, *m_dxt, *m_ddxr, *m_ddxt;
			double m_p_old, m_p_new, *m_p_data, m_dp;
			double *m_x_old, *m_x_new, *m_x_data, *m_dx;
			double *m_v_old, *m_v_new, *m_v_data, *m_dv;
			double *m_a_old, *m_a_new, *m_a_data, *m_da;
		};
	}
}