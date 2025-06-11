#pragma once

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

			void** m_args;
			convergence m_convergence;
			continuation m_continuation;
			stop_criteria m_stop_criteria;

			bool (*m_stop)(void**);
			void (*m_record)(void**);
			void (*m_update)(void**);
			void (*m_restore)(void**);
			void (*m_interface)(uint32_t, void**);

			uint32_t m_size, m_watch_dof;
			uint32_t m_step, m_attempt, m_iteration;
			uint32_t m_step_max, m_attempt_max, m_iteration_max;

			double *m_K, *m_C, *m_M;
			double *m_r, *m_fi, *m_fe;
			double m_p_old, m_p_new, *m_p_data;
			double *m_dxr, *m_dxt, *m_ddxr, *m_ddxt;
			double *m_x_old, *m_x_new, *m_x_data, m_dx;
			double *m_v_old, *m_v_new, *m_v_data, m_dv;
			double *m_a_old, *m_a_new, *m_a_data, m_da;
		};
	}
}