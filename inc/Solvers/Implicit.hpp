#pragma once

//Math
#include "Math/inc/Solvers/Solver.hpp"
#include "Math/inc/Solvers/Convergence.hpp"
#include "Math/inc/Solvers/Continuation.hpp"

namespace math
{
	namespace solvers
	{
		class Implicit : public virtual Solver
		{
		public:
			//constructor
			Implicit(void);

			//destructor
			~Implicit(void);

		protected:
			//solve
			void print(void) override;
			void setup(void) override;

			//compute
			virtual bool compute_equilibrium(void);

		public:
			//data
			bool m_equilibrium;
			Convergence m_convergence;
			Continuation m_continuation;
			uint32_t m_attempt, m_attempt_max;
			uint32_t m_iteration, m_iteration_max;
		};
	}
}