#pragma once

//std
#include <cstdint>
#include <initializer_list>

//math
#include "Math/Math/inc/linear/matrix.hpp"

namespace math
{
	class vector : public matrix
	{
	public:
		//constructors
		vector(void);
		vector(const matrix&);
		vector(const double*, uint32_t);
		vector(uint32_t, mode = mode::null);
		vector(std::initializer_list<double>);
		vector(double*, uint32_t, mode = mode::null);

		//destructor
		virtual ~vector(void);

		//size
		vector& resize(uint32_t);
		vector& resize(uint32_t, uint32_t);

		//linear
		vector& normalize(void);
		vector unit(void) const;
		matrix outer(void) const;
		matrix outer(const vector&) const;
		double inner(const vector&) const;
		double inner(const double*) const;
		

		//operators
		vector operator+(void) const;
		vector operator-(void) const;
		vector operator/(double) const;

		vector operator+(const vector&) const;
		vector operator-(const vector&) const;

		vector& operator=(double);
		vector& operator=(const double*);
		vector& operator=(const vector&);
		vector& operator=(std::initializer_list<double>);

		vector& operator+=(double);
		vector& operator-=(double);
		vector& operator*=(double);
		vector& operator/=(double);

		vector& operator+=(const double*);
		vector& operator-=(const double*);
		vector& operator+=(const vector&);
		vector& operator-=(const vector&);
	};
}