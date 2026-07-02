#pragma once

//Math
#include "Math/inc/Solvers/Solver.hpp"
#include "Math/inc/Solvers/StopCriteria.hpp"

namespace math
{
	namespace solvers
	{
		class StopCriteria;
	}
}

namespace math
{
	namespace solvers
	{
		class Incremental : public virtual Solver
		{
		public:
			//constructor
			Incremental(void);

			//destructor
			~Incremental(void);

			//solve
			void step(void);
			void solve(void) override;

			//serialization
			void save(const char*) const override;

		protected:
			//solve
			bool stop(void) override;
			void print(void) override;
			void setup(void) override;
			void record(void) override;
			
			//allocate
			void cleanup(void) override;
			void allocate_state(void) override;

		public:
			//data
			uint32_t m_step, m_step_max;
			StopCriteria m_stop_criteria;
			double *m_x_data, *m_v_data, *m_a_data, *m_p_data, *m_t_data;

			//friends
			friend class StopCriteria;
		};
	}
}