#pragma once

//std
#include <cstdint>
#include <functional>

//math
#include "Math/Math/inc/linear/matrix.hpp"
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/continuation.hpp"
#include "Math/Math/inc/solvers/stop_criteria.hpp"

namespace math
{
	namespace solvers
	{
		class newton_raphson
		{
		public:
			//constructors
			newton_raphson(void);

			//destructor
			~newton_raphson(void);

			//data
			void size(uint32_t);
			void save(const char*) const;

		private:
			//solve
			bool stop(void);
			void apply(void);
			void print(void);
			void setup(void);
			void update(void);
			void record(void);
			void restore(void);
			void predictor(void);
			void corrector(void);
			bool equilibrium(void);
			void load_predictor(void);
			void load_corrector(void);

		public:
			//solve
			void step(void);
			void solve(void);
			void cleanup(void);
			void allocate(void);

			//data
			bool m_silent;
			bool m_equilibrium;
			continuation m_continuation;
			stop_criteria m_stop_criteria;
			std::function<bool(void)> m_stop;
			std::function<void(void)> m_record;
			std::function<void(void)> m_update;
			std::function<void(void)> m_restore;
			std::function<void(uint32_t)> m_interface;
			std::function<void(double*, double, const double*)> m_residue;
			std::function<void(double*, double, const double*)> m_tangent_1;
			std::function<void(double*, double, const double*)> m_tangent_2;

			uint32_t m_watch_dof, m_size;
			uint32_t m_step, m_attempt, m_iteration;
			uint32_t m_step_max, m_attempt_max, m_iteration_max;
			double m_dp0, m_tolerance, m_p_new, *m_x_new, m_p_min, m_p_max, m_x_min, m_x_max;

		private:
			//data
			double m_p_old, *m_p_data, *m_x_old, *m_x_data;
			double *m_r, *m_g, *m_K, m_dp, m_ddp, *m_dx, *m_dxr, *m_dxt, *m_ddxr, *m_ddxt;

			//friends
			friend class stop_criteria;
		};
	}
}