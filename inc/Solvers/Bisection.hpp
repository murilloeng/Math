#pragma once

//std
#include <cstdint>

namespace math
{
	namespace solvers
	{
		class Bisection
		{
		public:
			//constructors
			Bisection(void);
	
			//destructor
			virtual ~Bisection(void);
	
		public:
			//solve
			bool solve(void**);
			bool solve(const void**);
	
			//data
			double m_x1, m_x2, m_xs, m_tolerance;
			uint32_t m_iteration, m_iteration_max;
	
			double(*m_system_1)(double, void**);
			double(*m_system_2)(double, const void**);
		};
	}
}