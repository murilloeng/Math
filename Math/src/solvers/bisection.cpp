//std
#include <cmath>
#include <cstdio>

//math
#include "Math/Math/inc/solvers/bisection.hpp"

namespace math
{
	//constructors
	bisection::bisection(void) : m_tolerance(1e-5), m_iteration_max(100), m_system_1(nullptr), m_system_2(nullptr)
	{
		return;
	}

	//destructor
	bisection::~bisection(void)
	{
		return;
	}
	
	//solve
	bool bisection::solve(void** args)
	{
		//data
		double fs;
		double f1 = m_system_1(m_x1, args);
		double f2 = m_system_1(m_x2, args);
		//check
		if(f1 * f2 > 0)
		{
			printf("Initial guess creates outputs with the same sign\n");
			return false;
		}
		//loop
		for(m_iteration = 0; m_iteration < m_iteration_max; m_iteration++)
		{
			//test
			if(fabs(f1) < m_tolerance) { m_xs = m_x1; return true; }
			if(fabs(f2) < m_tolerance) { m_xs = m_x2; return true; }
			//update
			m_xs = (m_x1 + m_x2) / 2;
			fs = m_system_1(m_xs, args);
			if(fs * f1 > 0) { m_x1 = m_xs; f1 = fs; } else { m_x2 = m_xs; f2 = fs; }
		}
		return false;
	}
	bool bisection::solve(const void** args)
	{
		//data
		double fs;
		double f1 = m_system_2(m_x1, args);
		double f2 = m_system_2(m_x2, args);
		//check
		if(f1 * f2 > 0)
		{
			printf("Initial guess creates outputs with the same sign\n");
			return false;
		}
		//loop
		for(m_iteration = 0; m_iteration < m_iteration_max; m_iteration++)
		{
			//test
			if(fabs(f1) < m_tolerance) { m_xs = m_x1; return true; }
			if(fabs(f2) < m_tolerance) { m_xs = m_x2; return true; }
			//update
			m_xs = (m_x1 + m_x2) / 2;
			fs = m_system_2(m_xs, args);
			if(fs * f1 > 0) { m_x1 = m_xs; f1 = fs; } else { m_x2 = m_xs; f2 = fs; }
		}
		return false;
	}
}