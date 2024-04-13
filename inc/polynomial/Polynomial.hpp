#pragma once

//std
#include <complex.h>

namespace math
{
	class Polynomial
	{
	public:
		//constructor
		Polynomial(void);

		//destructor
		~Polynomial(void);

		//data
		unsigned order(unsigned);
		unsigned order(void) const;

		double roots_error(void) const;
		const std::complex<double> root(unsigned) const;

		std::complex<double> coefficient(unsigned) const;
		std::complex<double> coefficient(unsigned, std::complex<double>);

		//value
		std::complex<double> value(std::complex<double>) const;
		std::complex<double> derivative(std::complex<double>) const;

		//analysis
		bool solve(void);

		//analysis
		void bound(void);
		void predictor(void);
		void corrector(void);

		//data
		double m_bound;
		unsigned m_order;
		double m_tolerance;
		unsigned m_iteration_max;
		std::complex<double>* m_roots;
		std::complex<double>* m_coefficients;
	};
}