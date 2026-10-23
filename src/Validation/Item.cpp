//std
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

//Math
#include "Math/inc/Validation/Item.hpp"

namespace math
{
	namespace validation
	{
		//constructor
		Item::Item(void) : m_silent{false}, m_tolerance{1e-2}
		{
			return;
		}

		//data
		void Item::tolerance(double tolerance)
		{
			m_tolerance = tolerance;
		}
		void Item::function(std::function<double(double)> function)
		{
			m_function = function;
		}
		void Item::adjuster_numeric(std::function<void(double&, double&)> adjuster_numeric)
		{
			m_adjuster_numeric = adjuster_numeric;
		}
		void Item::adjuster_reference(std::function<void(double&, double&)> adjuster_reference)
		{
			m_adjuster_reference = adjuster_reference;
		}

		//print
		void Item::print_numeric(void) const
		{
			for(const Point& point : m_points_numeric)
			{
				printf("%+.6e %+.6e\n", point.m_data[0], point.m_data[1]);
			}
		}
		void Item::print_reference(void) const
		{
			for(const Point& point : m_points_reference)
			{
				printf("%+.6e %+.6e\n", point.m_data[0], point.m_data[1]);
			}
		}

		//load
		void Item::load_numeric(const char* path, uint32_t col_1, uint32_t col_2)
		{
			load(m_points_numeric, path, col_1, col_2);
			if(m_adjuster_numeric) for(Point& point : m_points_numeric) m_adjuster_numeric(point.m_data[0], point.m_data[1]);
		}
		void Item::load_numeric(const double* x1, const double* x2, uint32_t rows)
		{
			load(m_points_numeric, x1, x2, rows);
			if(m_adjuster_numeric) for(Point& point : m_points_numeric) m_adjuster_numeric(point.m_data[0], point.m_data[1]);
		}
		void Item::load_numeric(const double* x, uint32_t rows, uint32_t cols, uint32_t col_1, uint32_t col_2)
		{
			load(m_points_numeric, x, rows, cols, col_1, col_2);
			if(m_adjuster_numeric) for(Point& point : m_points_numeric) m_adjuster_numeric(point.m_data[0], point.m_data[1]);
		}

		void Item::load_reference(const char* path, uint32_t col_1, uint32_t col_2)
		{
			load(m_points_reference, path, col_1, col_2);
			if(m_adjuster_reference) for(Point& point : m_points_reference) m_adjuster_reference(point.m_data[0], point.m_data[1]);
		}
		void Item::load_reference(const double* x1, const double* x2, uint32_t rows)
		{
			load(m_points_reference, x1, x2, rows);
			if(m_adjuster_reference) for(Point& point : m_points_reference) m_adjuster_reference(point.m_data[0], point.m_data[1]);
		}
		void Item::load_reference(const double* x, uint32_t rows, uint32_t cols, uint32_t col_1, uint32_t col_2)
		{
			load(m_points_reference, x, rows, cols, col_1, col_2);
			if(m_adjuster_reference) for(Point& point : m_points_reference) m_adjuster_reference(point.m_data[0], point.m_data[1]);
		}

		//bounds
		void Item::compute_bounds(void)
		{
			//setup
			m_bounds[0] = m_bounds[2] = +DBL_MAX;
			m_bounds[1] = m_bounds[3] = -DBL_MAX;
			//numeric
			for(const Point& point : m_points_numeric)
			{
				m_bounds[0] = fmin(m_bounds[0], point.m_data[0]);
				m_bounds[1] = fmax(m_bounds[1], point.m_data[0]);
				m_bounds[2] = fmin(m_bounds[2], point.m_data[1]);
				m_bounds[3] = fmax(m_bounds[3], point.m_data[1]);
			}
			//reference
			for(const Point& point : m_points_reference)
			{
				m_bounds[0] = fmin(m_bounds[0], point.m_data[0]);
				m_bounds[1] = fmax(m_bounds[1], point.m_data[0]);
				m_bounds[2] = fmin(m_bounds[2], point.m_data[1]);
				m_bounds[3] = fmax(m_bounds[3], point.m_data[1]);
			}
			//restore
			if(m_bounds[1] == m_bounds[0]) m_bounds[1] = m_bounds[0] + 1;
			if(m_bounds[3] == m_bounds[2]) m_bounds[3] = m_bounds[2] + 1;
		}
		Point Item::transform(const Point& point) const
		{
			return {
				(point.m_data[0] - m_bounds[0]) / (m_bounds[1] - m_bounds[0]),
				(point.m_data[1] - m_bounds[2]) / (m_bounds[3] - m_bounds[2])
			};
		}

		//load
		void Item::load(std::vector<Point>& points, const char* path, uint32_t col_1, uint32_t col_2)
		{
			//open
			char line[4096];
			FILE* file = fopen(path, "r");
			const uint32_t col_max = col_1 > col_2 ? col_1 : col_2;
			if(!file) throw std::runtime_error("Error: Item file " + std::string(path) + " could not be opened!");
			//read
			points.clear();
			while(fgets(line, sizeof(line), file))
			{
				Point point;
				uint32_t col = 0;
				char* end = nullptr;
				const char* ptr = line;
				for(col = 0; col <= col_max; col++)
				{
					const double value = strtod(ptr, &end);
					if(end == ptr) break;
					if(col == col_1) point.m_data[0] = value;
					if(col == col_2) point.m_data[1] = value;
					ptr = end;
				}
				if(col == 0) continue;
				if(col <= col_max)
				{
					fclose(file);
					throw std::runtime_error("Error: Item file " + std::string(path) + " has too few columns!");
				}
				points.push_back(point);
			}
			//close
			fclose(file);
		}
		void Item::load(std::vector<Point>& points, const double* x1, const double* x2, uint32_t rows)
		{
			points.resize(rows);
			for(uint32_t i = 0; i < rows; i++)
			{
				points[i].m_data[0] = x1[i];
				points[i].m_data[1] = x2[i];
			}
		}
		void Item::load(std::vector<Point>& points, const double* x, uint32_t rows, uint32_t cols, uint32_t col_1, uint32_t col_2)
		{
			points.resize(rows);
			for(uint32_t i = 0; i < rows; i++)
			{
				points[i].m_data[0] = x[i * cols + col_1];
				points[i].m_data[1] = x[i * cols + col_2];
			}
		}

		//validation
		bool Item::validate(void)
		{
			//check
			if(m_points_numeric.empty() || (!m_function && m_points_reference.empty()))
			{
				if(!m_silent) printf("Item has no data loaded!\n");
				return false;
			}
			//validate
			compute_bounds();
			return !m_function ? validate_data() : validate_function();
		}
		bool Item::validate_data(void)
		{
			//transform
			std::vector<Point> pn;
			pn.reserve(m_points_numeric.size());
			for(const Point& point : m_points_numeric) pn.push_back(transform(point));
			//validation
			m_error = 0;
			const uint64_t nr = m_points_reference.size();
			for(const Point& point_reference : m_points_reference)
			{
				const Point pr = transform(point_reference);
				double distance = pr.distance(pn[0]);
				for(uint32_t i = 0; i + 1 < pn.size(); i++)
				{
					distance = fmin(distance, pr.distance(pn[i], pn[i + 1]));
				}
				m_error += distance / nr;
			}
			if(m_error > m_tolerance && !m_silent)
			{
				printf("MAE: %+.2e\n", m_error);
			}
			//return
			return m_error < m_tolerance;
		}
		bool Item::validate_function(void)
		{
			//validation
			m_error = 0;
			const uint64_t nr = m_points_numeric.size();
			for(uint32_t i = 0; i < m_points_numeric.size(); i++)
			{
				const double x1 = m_points_numeric[i].m_data[0];
				const double x2 = m_points_numeric[i].m_data[1];
				m_error += fabs(x2 - m_function(x1)) / (m_bounds[3] - m_bounds[2]) / nr;
			}
			if(m_error > m_tolerance && !m_silent)
			{
				printf("MAE: %+.2e\n", m_error);
			}
			//return
			return m_error < m_tolerance;
		}
	}
}