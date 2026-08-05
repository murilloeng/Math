#pragma once

//std
#include <string>
#include <vector>
#include <cstdint>
#include <functional>

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
			void function(std::function<double(double)>);
			void adjuster_numeric(std::function<void(double&, double&)>);
			void adjuster_reference(std::function<void(double&, double&)>);

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
			//validation
			bool validate(void);
			bool validate_data(void);
			bool validate_function(void);

			//bounds
			void compute_bounds(void);
			Point transform(const Point&) const;

			//load
			void load(std::vector<Point>&, const char*, uint32_t, uint32_t);
			void load(std::vector<Point>&, const double*, const double*, uint32_t);
			void load(std::vector<Point>&, const double*, uint32_t, uint32_t, uint32_t, uint32_t);

			//data
			bool m_silent;
			double m_error;
			double m_tolerance;
			double m_bounds[4];
			std::vector<Point> m_points_numeric;
			std::vector<Point> m_points_reference;
			std::function<double(double)> m_function;
			std::function<void(double&, double&)> m_adjuster_numeric;
			std::function<void(double&, double&)> m_adjuster_reference;
		};
	}
}