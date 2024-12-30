#pragma once

//std
#include <cstdint>
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
		uint32_t order(uint32_t);
		uint32_t order(void) const;

		double roots_error(void) const;
		const std::complex<double> root(uint32_t) const;

		std::complex<double> coefficient(uint32_t) const;
		std::complex<double> coefficient(uint32_t, std::complex<double>);

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
		uint32_t m_order;
		double m_tolerance;
		uint32_t m_iteration_max;
		std::complex<double>* m_roots;
		std::complex<double>* m_coefficients;
	};
}