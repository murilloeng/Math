//std
#include <cmath>

//Math
#include "Math/inc/Validation/Point.hpp"

namespace math
{
	namespace validation
	{
		//constructor
		Point::Point(double x1, double x2) : m_data{x1, x2}
		{
			return;
		}

		//destructor
		Point::~Point(void)
		{
			return;
		}

		//algebra
		double Point::norm(void) const
		{
			return sqrt(inner(*this));
		}
		double Point::inner(const Point& p) const
		{
			return m_data[0] * p.m_data[0] + m_data[1] * p.m_data[1];
		}

		//operators
		Point Point::operator+(const Point& p) const
		{
			return {m_data[0] + p.m_data[0], m_data[1] + p.m_data[1]};
		}
		Point Point::operator-(const Point& p) const
		{
			return {m_data[0] - p.m_data[0], m_data[1] - p.m_data[1]};
		}
		Point operator*(double s, const Point& p)
		{
			return {s * p.m_data[0], s * p.m_data[1]};
		}

		//distance
		double Point::distance(const Point& p) const
		{
			return (*this - p).norm();
		}
		double Point::distance(const Point& p1, const Point& p2) const
		{
			const Point pr = p2 - p1;
			const double t = (*this - p1).inner(pr) / pr.inner(pr);
			return (p1 + fmax(fmin(t, 1), 0) * pr - *this).norm();
		}
	}
}