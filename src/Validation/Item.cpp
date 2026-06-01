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
		Item::Item(void) : m_tolerance{1e-5}, m_function{nullptr}
		{
			return;
		}

		//destructor
		Item::~Item(void)
		{
			return;
		}

		//distance
		double Item::distance(void) const
		{
			return m_distance;
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
			return !m_function ? validate_data() : validate_function();
		}
		bool Item::validate_data(void)
		{
			//validation
			for(uint32_t i = 0; i < m_points_reference.size(); i++)
			{
				bool test = false;
				m_distance = DBL_MAX;
				for(uint32_t j = 0; j + 1 < m_points_numeric.size(); j++)
				{
					if(!test)
					{
						m_distance = fmin(m_distance, m_points_reference[i].distance(m_points_numeric[j], m_points_numeric[j + 1]));
						test = m_distance < m_tolerance;
					}
				}
				if(!test)
				{
					printf("minimal distance: %+.2e\n", m_distance);
					return false;
				}
			}
			//return
			return true;
		}
		bool Item::validate_function(void)
		{
			//validation
			bool test = true;
			for(uint32_t i = 0; i < m_points_numeric.size(); i++)
			{
				const double x1 = m_points_numeric[i].m_data[0];
				const double x2 = m_points_numeric[i].m_data[1];
				test = test && fabs(x2 - m_function(x1)) < m_tolerance;
			}
			//return
			return test;
		}
	}
}