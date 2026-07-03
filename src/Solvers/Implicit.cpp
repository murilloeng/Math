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