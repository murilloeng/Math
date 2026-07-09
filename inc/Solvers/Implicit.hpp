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
		protected:
			//constructor
			Implicit(void);

			//destructor
			~Implicit(void);

		public:
			//data
			uint32_t attempt(void) const;
			uint32_t iteration(void) const;

			uint32_t attempt_max(uint32_t);
			uint32_t attempt_max(void) const;

			uint32_t iteration_max(uint32_t);
			uint32_t iteration_max(void) const;

			Convergence& convergence(void);
			Continuation& continuation(void);
		
		protected:
			//solve
			void print(void) override;
			void setup(void) override;

			//compute
			virtual bool equilibrium(void);

			//data
			Convergence m_convergence;
			Continuation m_continuation;
			uint32_t m_attempt, m_attempt_max;
			uint32_t m_iteration, m_iteration_max;
		};
	}
}