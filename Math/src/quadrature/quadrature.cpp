//std
#include <cstdlib>

//math
#include "Math/Math/inc/Quadrature/Rule.hpp"
#include "Math/Math/inc/Quadrature/Quadrature.hpp"

extern "C"
{
	void lobatto_set(int32_t, double[], double[]);
	void legendre_set(int32_t, double[], double[]);
}

namespace math
{
	namespace quadrature
	{
		//constructor
		Quadrature::Quadrature(uint32_t order) : m_points(nullptr), m_weights(nullptr), m_rule(rule::legendre)
		{
			this->order(order);
		}
		Quadrature::Quadrature(quadrature::rule rule, uint32_t order) : m_points(nullptr), m_weights(nullptr), m_rule(rule)
		{
			this->order(order);
		}

		//destructor
		Quadrature::~Quadrature(void)
		{
			delete[] m_points;
			delete[] m_weights;
		}

		//serialization
		void Quadrature::load(FILE* file)
		{
			uint32_t rule;
			if(fscanf(file, "%d %d", &rule, &m_order) != 2)
			{
				printf("\tError: Unable to load Quadrature!\n");
				exit(EXIT_FAILURE);
			}
			m_rule = quadrature::rule(rule);
		}
		void Quadrature::save(FILE* file) const
		{
			fprintf(file, "%02d %02d ", (uint32_t) m_rule, m_order);
		}

		//data
		uint32_t Quadrature::order(void) const
		{
			return m_order;
		}
		uint32_t Quadrature::order(uint32_t order)
		{
			m_order = order;
			delete[] m_points;
			delete[] m_weights;
			m_points = new double[m_order];
			m_weights = new double[m_order];
			if(m_rule == quadrature::rule::lobatto)
			{
				lobatto_set(m_order, m_points, m_weights);
			}
			else if(m_rule == quadrature::rule::legendre)
			{
				legendre_set(m_order, m_points, m_weights);
			}
			return m_order;
		}

		quadrature::rule Quadrature::rule(void) const
		{
			return m_rule;
		}
		quadrature::rule Quadrature::rule(quadrature::rule rule)
		{
			return m_rule = rule;
		}

		//name
		const char* Quadrature::rule_name(void) const
		{
			return rule_name(m_rule);
		}
		const char* Quadrature::rule_name(quadrature::rule rule)
		{
			return 
				rule == quadrature::rule::lobatto ? "Lobatto" :
				rule == quadrature::rule::legendre ? "Legendre" : "Error";
		}

		//points
		const double* Quadrature::points(void) const
		{
			return m_points;
		}
		const double* Quadrature::weights(void) const
		{
			return m_weights;
		}
		double Quadrature::point(uint32_t index) const
		{
			return m_points[index];
		}
		double Quadrature::weight(uint32_t index) const
		{
			return m_weights[index];
		}
	}
}