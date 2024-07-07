#pragma once

namespace math
{
	class bisection
	{
	public:
		//constructors
		bisection(void);

		//destructor
		virtual ~bisection(void);

	public:
		//solve
		bool solve(void**);

		//data
		double(*m_system)(double, void**);
		double m_x1, m_x2, m_xs, m_tolerance;
		unsigned m_iteration, m_iteration_max;
	};
}