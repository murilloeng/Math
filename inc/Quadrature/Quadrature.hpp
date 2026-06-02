#pragma once

//std
#include <cstdio>
#include <cstdint>

namespace math
{
	namespace quadrature
	{
		enum class Rule : uint32_t;
	}
}

namespace math
{
	namespace quadrature
	{
		class Quadrature
		{
		public:
			//constructor
			Quadrature(uint32_t);
			Quadrature(quadrature::Rule, uint32_t);

			//destructor
			~Quadrature(void);

			//serialization
			void load(FILE*);
			void save(FILE*) const;

			//data
			uint32_t order(uint32_t);
			uint32_t order(void) const;

			quadrature::Rule rule(void) const;
			quadrature::Rule rule(quadrature::Rule);

			//name
			const char* rule_name(void) const;
			static const char* rule_name(quadrature::Rule);

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
			quadrature::Rule m_rule;
		};
	}
}