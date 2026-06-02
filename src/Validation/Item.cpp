//std
#include <cmath>
#include <cfloat>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Validation/Item.hpp"

namespace math
{
	namespace validation
	{
		//constructor
		Item::Item(void) : m_silent{false}, m_tolerance{1e-2}, m_function{nullptr}
		{
			return;
		}

		//destructor
		Item::~Item(void)
		{
			return;
		}

		//data
		void Item::tolerance(double tolerance)
		{
			m_tolerance = tolerance;
		}
		void Item::function(double(*function)(double))
		{
			m_function = function;
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
		}
		void Item::load_numeric(const double* x1, const double* x2, uint32_t rows)
		{
			load(m_points_numeric, x1, x2, rows);
		}
		void Item::load_numeric(const double* x, uint32_t rows, uint32_t cols, uint32_t col_1, uint32_t col_2)
		{
			load(m_points_numeric, x, rows, cols, col_1, col_2);
		}

		void Item::load_reference(const char* path, uint32_t col_1, uint32_t col_2)
		{
			load(m_points_reference, path, col_1, col_2);
		}
		void Item::load_reference(const double* x1, const double* x2, uint32_t rows)
		{
			load(m_points_reference, x1, x2, rows);
		}
		void Item::load_reference(const double* x, uint32_t rows, uint32_t cols, uint32_t col_1, uint32_t col_2)
		{
			load(m_points_reference, x, rows, cols, col_1, col_2);
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
			char line[512];
			uint32_t rows = 0, cols = 0;
			FILE* file = fopen(path, "r");
			//rows
			while(!feof(file))
			{
				if(fgetc(file) == '\n') rows++;
			}
			//cols
			rewind(file);
			if(!fgets(line, sizeof(line), file))
			{
				throw std::runtime_error("Error: Item file read failed!");
			}
			for(char c : line)
			{
				if(c == ' ') cols++; else if(c == '\n') break;
			}
			//read
			rewind(file);
			points.resize(rows);
			for(uint32_t i = 0; i < rows; i++)
			{
				for(uint32_t j = 0; j < cols; j++)
				{
					if(j != col_1 && j != col_2)
					{
						if(fscanf(file, "%*f") != 0) throw std::runtime_error("Error: Item file read failed!");
					}
					else if(j == col_1)
					{
						if(fscanf(file, "%lf", &points[i].m_data[0]) != 1) throw std::runtime_error("Error: Item file read failed!");
					}
					else if(j == col_2)
					{
						if(fscanf(file, "%lf", &points[i].m_data[1]) != 1) throw std::runtime_error("Error: Item file read failed!");
					}
				}
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
			compute_bounds();
			return !m_function ? validate_data() : validate_function();
		}
		bool Item::validate_data(void)
		{
			//validation
			m_error = 0;
			const uint64_t nr = m_points_reference.size();
			for(const Point& point_reference : m_points_reference)
			{
				double distance = DBL_MAX;
				const Point pr = transform(point_reference);
				for(const Point& point_numeric : m_points_numeric)
				{
					const Point& pn = transform(point_numeric);
					distance = fmin(distance, pr.distance(pn));
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