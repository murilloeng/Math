//std
#include <cstdlib>

//math
#include "Math/Math/inc/quadrature/rule.hpp"
#include "Math/Math/inc/quadrature/quadrature.hpp"

extern "C"
{
	void lobatto_set(int, double[], double[]);
	void legendre_set(int, double[], double[]);
}

namespace math
{
	//constructor
	quadrature::quadrature(uint32_t order) : m_points(nullptr), m_weights(nullptr), m_rule(rule::legendre)
	{
		this->order(order);
	}
	quadrature::quadrature(math::rule rule, uint32_t order) : m_points(nullptr), m_weights(nullptr), m_rule(rule)
	{
		this->order(order);
	}

	//destructor
	quadrature::~quadrature(void)
	{
		delete[] m_points;
		delete[] m_weights;
	}

	//serialization
	void quadrature::load(FILE* file)
	{
		uint32_t rule;
		if(fscanf(file, "%d %d", &rule, &m_order) != 2)
		{
			printf("\tError: Unable to load Quadrature!\n");
			exit(EXIT_FAILURE);
		}
		m_rule = math::rule(rule);
	}
	void quadrature::save(FILE* file) const
	{
		fprintf(file, "%02d %02d ", (uint32_t) m_rule, m_order);
	}

	//data
	uint32_t quadrature::order(void) const
	{
		return m_order;
	}
	uint32_t quadrature::order(uint32_t order)
	{
		m_order = order;
		delete[] m_points;
		delete[] m_weights;
		m_points = new double[m_order];
		m_weights = new double[m_order];
		if(m_rule == math::rule::lobatto)
		{
			lobatto_set(m_order, m_points, m_weights);
		}
		else if(m_rule == math::rule::legendre)
		{
			legendre_set(m_order, m_points, m_weights);
		}
		return m_order;
	}

	math::rule quadrature::rule(void) const
	{
		return m_rule;
	}
	math::rule quadrature::rule(math::rule rule)
	{
		return m_rule = rule;
	}

	//name
	const char* quadrature::rule_name(void) const
	{
		return rule_name(m_rule);
	}
	const char* quadrature::rule_name(math::rule rule)
	{
		return 
			rule == math::rule::lobatto ? "Lobatto" :
			rule == math::rule::legendre ? "Legendre" : "Error";
	}

	//points
	const double* quadrature::points(void) const
	{
		return m_points;
	}
	const double* quadrature::weights(void) const
	{
		return m_weights;
	}
	double quadrature::point(uint32_t index) const
	{
		return m_points[index];
	}
	double quadrature::weight(uint32_t index) const
	{
		return m_weights[index];
	}
}