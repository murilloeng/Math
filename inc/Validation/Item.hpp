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
			//data
			void tolerance(double);
			void function(double(*)(double));

			void data_numeric(FILE*);
			void data_numeric(FILE*, uint32_t);
			void data_numeric(const double*, uint32_t);
			void data_numeric(const double*, uint32_t, uint32_t, uint32_t);

			void data_reference(FILE*);
			void data_reference(FILE*, uint32_t);
			void data_reference(const double*, uint32_t);
			void data_reference(const double*, uint32_t, uint32_t, uint32_t);
	
			//distance
			double distance(void) const;
	
		private:
			//data
			void data(std::vector<Point>&, const char*, uint32_t, uint32_t);
			void data(std::vector<Point>&, const double*, const double*, uint32_t);
			void data(std::vector<Point>&, const double*, uint32_t, uint32_t, uint32_t, uint32_t);

			//validation
			bool validate(void);
	
			//load
			void load_numeric(void);
			void load_reference(FILE*);
	
			//validation
			bool validate_file(void);
			bool validate_function(void);
	
			//data
			double m_distance;
			double m_tolerance;
			double(*m_function)(double);
			std::vector<Point> m_data_numeric;
			std::vector<Point> m_data_reference;
		};
	}
}