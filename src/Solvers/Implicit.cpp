//std
#include <cstdio>

//Math
#include "Math/inc/Solvers/Implicit.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		Implicit::Implicit(void) :
			m_convergence{this}, m_continuation{this},
			m_attempt{0}, m_attempt_max{5}, m_iteration{0}, m_iteration_max{10}
		{
			return;
		}
		
		//destructor
		Implicit::~Implicit(void)
		{
			return;
		}

		//data
		uint32_t Implicit::attempt(void) const
		{
			return m_attempt;
		}
		uint32_t Implicit::iteration(void) const
		{
			return m_iteration;
		}

		uint32_t Implicit::attempt_max(void) const
		{
			return m_attempt_max;
		}
		uint32_t Implicit::attempt_max(uint32_t attempt_max)
		{
			return m_attempt_max;
		}

		uint32_t Implicit::iteration_max(void) const
		{
			return m_iteration_max;
		}
		uint32_t Implicit::iteration_max(uint32_t iteration_max)
		{
			return m_iteration_max;
		}

		Convergence& Implicit::convergence(void)
		{
			return m_convergence;
		}
		Continuation& Implicit::continuation(void)
		{
			return m_continuation;
		}

		//solve
		void Implicit::print(void)
		{
			if(m_silent) return;
			printf("Attempts: %2d ", m_attempt);
			printf("Iterations: %2d ", m_iteration);
		}
		void Implicit::setup(void)
		{
			m_attempt = 0;
			m_iteration = 0;
		}

		//analysis
		bool Implicit::equilibrium(void)
		{
			return m_status = m_convergence.check();
		}
	}
}