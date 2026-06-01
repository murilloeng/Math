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

		//data
		void Item::tolerance(double tolerance)
		{
			m_tolerance = tolerance;
		}
		void Item::function(double(*function)(double))
		{
			m_function = function;
		}

		void Item::data_numeric(FILE*)
		{

		}
		void Item::data_numeric(FILE*, uint32_t)

		{
			
		}
		void Item::data_numeric(const double*, uint32_t)

		{
			
		}
		void Item::data_numeric(const double*, uint32_t, uint32_t, uint32_t)
		{
			
		}

		void Item::data_reference(FILE*)
		{
			
		}
		void Item::data_reference(FILE*, uint32_t)
		{
			
		}
		void Item::data_reference(const double*, uint32_t)
		{
			
		}
		void Item::data_reference(const double*, uint32_t, uint32_t, uint32_t)
		{
			
		}

		//distance
		double Item::distance(void) const
		{
			return m_distance;
		}

		//data
		void Item::data(std::vector<Point>& data, const char* path, uint32_t index_1, uint32_t index_2)
		{
			//open
			char line[2048];
			uint32_t rows = 0, cols = 0;
			FILE* file = fopen(path, "r");
			//rows
			while(!feof(file))
			{
				if(fgetc(file) == '\n') rows++;
			}
			//cols
			rewind(file);
			fgets(line, sizeof(line), file);
			for(char c : line)
			{
				if(c == ' ') cols++; else if(c == '\n') break;
			}
			//read
			rewind(file);
			data.resize(rows);
			for(uint32_t i = 0; i < rows; i++)
			{
				for(uint32_t j = 0; j < cols; j++)
				{
					if(j != index_1 && j != index_2) fscanf(file, "%*lf");
					else if(j == index_1) fscanf(file, "%lf", &data[i].m_data[0]);
					else if(j == index_2) fscanf(file, "%lf", &data[i].m_data[1]);
				}
			}
			//close
			fclose(file);
		}
		void Item::data(std::vector<Point>& data, const double* x1, const double* x2, uint32_t rows)
		{
			data.resize(rows);
			for(uint32_t i = 0; i < rows; i++)
			{
				data[i].m_data[0] = x1[i];
				data[i].m_data[1] = x2[i];
			}
		}
		void Item::data(std::vector<Point>& data, const double*, uint32_t, uint32_t, uint32_t, uint32_t)
		{

		}

		//validation
		bool Item::validate(void)
		{
			load_numeric();
			return !m_function ? validate_file() : validate_function();
		}

		//load
		void Item::load_numeric(void)
		{
			// //data
			// const uint32_t ss = m_model->analysis()->solver()->state_set();
			// const uint32_t steps = m_model->analysis()->solver()->step() + 1;
			// const uint32_t ps = ss & uint32_t(fea::analysis::solvers::state::p);
			// //load
			// m_data_numeric.resize(steps);
			// for(uint32_t i = 0; i < steps; i++)
			// {
			// 	m_data_numeric[i].m_data[0] = 
			// 		m_nodes[1] == UINT32_MAX && !ps ? 
			// 		m_model->analysis()->solver()->state_data(i) :
			// 		m_model->mesh()->node(m_nodes[0])->state_data(m_dof[0], i);
			// 	m_data_numeric[i].m_data[1] = 
			// 		m_nodes[1] == UINT32_MAX ? 
			// 		ps ? m_model->analysis()->solver()->state_data(i) : 
			// 		m_model->mesh()->node(m_nodes[0])->state_data(m_dof[0], i) :
			// 		m_model->mesh()->node(m_nodes[1])->state_data(m_dof[1], i) ;
			// }
		}
		void Item::load_reference(FILE* file)
		{
			//data
			uint32_t lines = 0;
			//count
			while(!feof(file))
			{
				if(fgetc(file) == '\n') lines++;
			}
			//load
			rewind(file);
			m_data_reference.resize(lines);
			for(size_t i = 0; i < lines; i++)
			{
				fscanf(file, "%lf %lf", &m_data_reference[i].m_data[0], &m_data_reference[i].m_data[1]);
			}
		}

		bool Item::validate_file(void)
		{
			// //data
			// char buffer[512];
			// if(!strlen(m_validator_path))
			// {
			// 	strcpy(buffer, m_path.c_str());
			// }
			// else
			// {
			// 	sprintf(buffer, "%s/%s", m_validator_path, m_path.c_str());
			// }
			// FILE* file = fopen(buffer, "r");
			// const uint32_t steps = m_model->analysis()->solver()->step();
			// //check
			// if(!file)
			// {
			// 	throw std::runtime_error("Unable to open validation file!");
			// }
			// //validation
			// load_reference(file);
			// for(uint32_t i = 0; i < m_data_reference.size(); i++)
			// {
			// 	bool test = false;
			// 	m_distance = DBL_MAX;
			// 	for(uint32_t j = 1; j < m_data_numeric.size(); j++)
			// 	{
			// 		if(!test)
			// 		{
			// 			m_distance = fmin(m_distance, m_data_reference[i].distance(m_data_numeric[j - 1], m_data_numeric[j]));
			// 			test = m_distance < m_tolerance;
			// 		}
			// 	}
			// 	if(!test)
			// 	{
			// 		return false;
			// 	}
			// }
			// fclose(file);
			//return
			return true;
		}
		bool Item::validate_function(void)
		{
			//validation
			bool test = true;
			for(uint32_t i = 0; i < m_data_numeric.size(); i++)
			{
				const double x1 = m_data_numeric[i].m_data[0];
				const double x2 = m_data_numeric[i].m_data[1];
				test = test && fabs(x2 - m_function(x1)) < m_tolerance;
			}
			//return
			return test;
		}
	}
}