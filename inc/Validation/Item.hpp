#pragma once

//std
#include <string>
#include <vector>
#include <cstdint>

//Math
#include "Math/inc/Validation/Point.hpp"

namespace math
{
	namespace validation
	{
		class Validator;
	}
}

namespace math
{
	namespace validation
	{
		class Item
		{
			//friends
			friend class Validator;
	
		private:
			//constructor
			Item(void);
	
			//destructor
			~Item(void);
	
		public:
			//distance
			double distance(void) const;

			//data
			void tolerance(double);
			void function(double(*)(double));

			//print
			void print_numeric(void) const;
			void print_reference(void) const;

			//load
			void load_numeric(const char*, uint32_t, uint32_t);
			void load_numeric(const double*, const double*, uint32_t);
			void load_numeric(const double*, uint32_t, uint32_t, uint32_t, uint32_t);

			void load_reference(const char*, uint32_t, uint32_t);
			void load_reference(const double*, const double*, uint32_t);
			void load_reference(const double*, uint32_t, uint32_t, uint32_t, uint32_t);
	
		private:
			//load
			void load(std::vector<Point>&, const char*, uint32_t, uint32_t);
			void load(std::vector<Point>&, const double*, const double*, uint32_t);
			void load(std::vector<Point>&, const double*, uint32_t, uint32_t, uint32_t, uint32_t);

			//validation
			bool validate(void);
			bool validate_data(void);
			bool validate_function(void);
	
			//data
			double m_distance;
			double m_tolerance;
			double(*m_function)(double);
			std::vector<Point> m_points_numeric;
			std::vector<Point> m_points_reference;
		};
	}
}