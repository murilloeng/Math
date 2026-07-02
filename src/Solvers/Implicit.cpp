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
			m_equilibrium{false}, m_convergence{this}, m_continuation{this, Continuation::Type::ArcLengthCylindrical},
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
		Convergence& Implicit::convergence(void)
		{
			return m_convergence;
		}
		Continuation& Implicit::continuation(void)
		{
			return m_continuation;
		}

		uint32_t Implicit::attempt_max(void) const
		{
			return m_attempt_max;
		}
		uint32_t Implicit::attempt_max(uint32_t attempt_max)
		{
			return m_attempt_max = attempt_max;
		}

		uint32_t Implicit::iteration_max(void) const
		{
			return m_iteration_max;
		}
		uint32_t Implicit::iteration_max(uint32_t iteration_max)
		{
			return m_iteration_max;
		}

		//analysis
		bool Implicit::equilibrium(void)
		{
			return m_equilibrium = m_convergence.check();
		}

		//solve
		void Implicit::print(void)
		{
			printf("Attempts: %2d ", m_attempt);
			printf("Iterations: %2d ", m_iteration);
		}
		void Implicit::setup(void)
		{
			m_attempt = 0;
			m_iteration = 0;
		}
	}
}