//std
#include <cmath>

//roots
#include "Math/inc/polynomial/Polynomial.hpp"

namespace math
{
	//constructor
	Polynomial::Polynomial(void) : m_bound(0), m_order(0), m_tolerance(1e-5), m_iteration_max(10000), m_roots(nullptr), m_coefficients(nullptr)
	{
		return;
	}

	//destructor
	Polynomial::~Polynomial(void)
	{
		delete[] m_roots;
		delete[] m_coefficients;
	}

	//data
	unsigned Polynomial::order(void) const
	{
		return m_order;
	}
	unsigned Polynomial::order(unsigned order)
	{
		delete[] m_roots;
		delete[] m_coefficients;
		m_roots = new std::complex<double>[order];
		m_coefficients = new std::complex<double>[order + 1];
		return m_order = order;
	}

	double Polynomial::roots_error(void) const
	{
		double error = 0;
		for(unsigned i = 0; i < m_order; i++)
		{
			std::complex<double> v = value(m_roots[i]);
			error = fmax(error, sqrt(v.real() * v.real() + v.imag() * v.imag()));
		}
		return error;
	}
	const std::complex<double> Polynomial::root(unsigned index) const
	{
		return m_roots[index];
	}

	std::complex<double> Polynomial::coefficient(unsigned index) const
	{
		return m_coefficients[index];
	}
	std::complex<double> Polynomial::coefficient(unsigned index, std::complex<double> coefficient)
	{
		return m_coefficients[index] = coefficient;
	}

	//value
	std::complex<double> Polynomial::value(std::complex<double> x) const
	{
		//data
		std::complex<double> v(0, 0), e(1, 0);
		//value
		for(unsigned i = 0; i <= m_order; i++)
		{
			if(i != 0) e *= x;
			v += m_coefficients[i] * e;
		}
		//return
		return v;
	}
	std::complex<double> Polynomial::derivative(std::complex<double> x) const
	{
		//data
		std::complex<double> v(0, 0), e(1, 0);
		//value
		for(unsigned i = 1; i <= m_order; i++)
		{
			if(i != 1) e *= x;
			v += double(i) * m_coefficients[i] * e;
		}
		//return
		return v;
	}

	//analysis
	bool Polynomial::solve(void)
	{
		bound();
		predictor();
		for(unsigned i = 0; i < m_iteration_max; i++)
		{
			corrector();
			if(roots_error() < m_tolerance)
			{
				return true;
			}
		}
		return false;
	}
	void Polynomial::bound(void)
	{
		double b1 = 0, b2 = 0;
		for(unsigned i = 0; i < m_order; i++)
		{
			const std::complex<double> c = m_coefficients[i] / m_coefficients[m_order];
			const double s = sqrt(c.real() * c.real() + c.imag() * c.imag());
			b1 += s;
			b2 = fmax(b2, s);
		}
		b2 += 1;
		b1 = fmax(1, b1);
		m_bound = fmin(b1, b2);
	}
	void Polynomial::predictor(void)
	{
		for(unsigned i = 0; i < m_order; i++)
		{
			const double r = m_bound * double(rand()) / RAND_MAX;
			const double t = 2 * M_PI * double(rand()) / RAND_MAX;
			m_roots[i] = std::complex<double>(r * cos(t), r * sin(t));
		}
	}
	void Polynomial::corrector(void)
	{
		for(unsigned i = 0; i < m_order; i++)
		{
			std::complex<double> s(0, 0);
			std::complex<double> v = value(m_roots[i]);
			std::complex<double> d = derivative(m_roots[i]);
			for(unsigned j = 0; j < m_order; j++)
			{
				if(i != j) s += 1.0 / (m_roots[i] - m_roots[j]);
			}
			m_roots[i] = m_roots[i] - v / (d - v * s);
		}
	}
}
