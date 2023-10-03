#pragma once

//std
#include <cstdio>

namespace math
{
	enum class rule : unsigned;
}

namespace math
{
	class quadrature
	{
	public:
		//constructor
		quadrature(unsigned);
		quadrature(math::rule, unsigned);

		//destructor
		~quadrature(void);

		//serialization
		void load(FILE*);
		void save(FILE*) const;

		//data
		unsigned order(unsigned);
		unsigned order(void) const;

		math::rule rule(void) const;
		math::rule rule(math::rule);

		//name
		const char* rule_name(void) const;
		static const char* rule_name(math::rule);

		//points
		double point(unsigned) const;
		double weight(unsigned) const;
		const double* points(void) const;
		const double* weights(void) const;

	private:
		//data
		unsigned m_order;
		double* m_points;
		double* m_weights;
		math::rule m_rule;
	};
}