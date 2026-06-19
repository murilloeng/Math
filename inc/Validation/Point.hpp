#pragma once

namespace math
{
	namespace validation
	{
		struct Point
		{
			//constructor
			Point(double = 0, double = 0);

			//destructor
			~Point(void);

			//algebra
			double norm(void) const;
			double inner(const Point&) const;

			//operators
			Point operator+(const Point&) const;
			Point operator-(const Point&) const;
			friend Point operator*(double, const Point&);

			//distance
			double distance(const Point&) const;
			double distance(const Point&, const Point&) const;

			//data
			double m_data[2];
		};
	}
}