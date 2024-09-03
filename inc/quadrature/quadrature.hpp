#pragma once

//std
#include <cstdio>
#include <cstdint>

namespace math
{
	enum class rule : uint32_t;
}

namespace math
{
	class quadrature
	{
	public:
		//constructor
		quadrature(uint32_t);
		quadrature(math::rule, uint32_t);

		//destructor
		~quadrature(void);

		//serialization
		void load(FILE*);
		void save(FILE*) const;

		//data
		uint32_t order(uint32_t);
		uint32_t order(void) const;

		math::rule rule(void) const;
		math::rule rule(math::rule);

		//name
		const char* rule_name(void) const;
		static const char* rule_name(math::rule);

		//points
		double point(uint32_t) const;
		double weight(uint32_t) const;
		const double* points(void) const;
		const double* weights(void) const;

	private:
		//data
		uint32_t m_order;
		double* m_points;
		double* m_weights;
		math::rule m_rule;
	};
}