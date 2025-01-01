//std
#include <cmath>
#include <cstdio>
#include <cstring>
#include <malloc.h>

//ext
#include "external/cpp/inc/suitesparse/umfpack.h"

//math
#include "Math/Math/inc/linear/sparse.hpp"
#include "Math/Math/inc/linear/vector.hpp"

namespace math
{
	//constructors
	sparse::sparse(double* ptr, const int32_t* row_map, const int32_t* col_map, uint32_t rows, uint32_t cols) : 
		m_rows(rows), m_cols(cols), m_ptr(ptr), m_ref(ptr), m_row_map(row_map), m_col_map(col_map)
	{
		return;
	}
	sparse::sparse(const double* ref, const int32_t* row_map, const int32_t* col_map, uint32_t rows, uint32_t cols) : 
		m_rows(rows), m_cols(cols), m_ptr(nullptr), m_ref(ref), m_row_map(row_map), m_col_map(col_map)
	{
		return;
	}

	//constructors
	sparse::~sparse(void)
	{
		return;
	}

	//data
	double* sparse::data(void)
	{
		return m_ptr;
	}
	const double* sparse::data(void) const
	{
		return m_ref;
	}

	//linear
	double sparse::norm(void) const
	{
		double s = 0;
		for(uint32_t i = 0; i < m_cols; i++)
		{
			for(int32_t j = m_col_map[i]; j < m_col_map[i + 1]; j++)
			{
				s += m_ref[j] * m_ref[j];
			}
		}
		return sqrt(s);
	}
	double sparse::trace(void) const
	{
		double s = 0;
		for(uint32_t i = 0; i < m_cols; i++)
		{
			for(int32_t j = m_col_map[i]; j < m_col_map[i + 1]; j++)
			{
				if(i == (uint32_t) m_row_map[j])
				{
					s += m_ref[j];
				}
			}
		}
		return s;
	}

	//print
	void sparse::print(const char* header, bool dense) const
	{
		//header
		if(strlen(header) != 0)
		{
			printf("%s\n", header);
		}
		//print
		if(!dense)
		{
			for(uint32_t i = 0; i < m_cols; i++)
			{
				for(int32_t j = m_col_map[i]; j < m_col_map[i + 1]; j++)
				{
					printf("(%04d, %04d): %+.2e\n", m_row_map[j], i, m_ref[j]);
				}
			}
		}
		else
		{
			for(uint32_t i = 0; i < m_rows; i++)
			{
				for(uint32_t j = 0; j < m_cols; j++)
				{
					printf("%+.2e ", (*this)(i, j));
				}
				printf("\n");
			}
		}
	}

	//operators
	vector sparse::operator*(const vector& v) const
	{
		//check
		if(m_cols != v.rows())
		{
			fprintf(stderr, "Error: Incompatible multiplication of sparse matrix and vector!\n");
			exit(EXIT_FAILURE);
		}
		//vector
		vector r(v.rows(), mode::zeros);
		for(uint32_t i = 0; i < m_cols; i++)
		{
			for(int32_t j = m_col_map[i]; j < m_col_map[i + 1]; j++)
			{
				r[m_row_map[j]] += v[i] * m_ref[j];
			}
		}
		return r;
	}
	double& sparse::operator()(uint32_t i, uint32_t j)
	{
		for(int32_t p = m_col_map[j]; p < m_col_map[j + 1]; p++)
		{
			if(int32_t(i) == m_row_map[p])
			{
				return m_ptr[p];
			}
		}
		return m_ptr[0];
	}
	const double& sparse::operator()(uint32_t i, uint32_t j) const
	{
		for(int32_t p = m_col_map[j]; p < m_col_map[j + 1]; p++)
		{
			if(int32_t(i) == m_row_map[p])
			{
				return m_ref[p];
			}
		}
		fprintf(stderr, "Error: Sparse matrix operator() called with out of range index!\n");
		exit(EXIT_FAILURE);
	}

	//convert
	matrix sparse::convert(void) const
	{
		matrix M(m_rows, m_cols, mode::zeros);
		for(uint32_t i = 0; i < m_cols; i++)
		{
			for(int32_t j = m_col_map[i]; j < m_col_map[i + 1]; j++)
			{
				M(m_row_map[j], i) = m_ref[j];
			}
		}
		return M;
	}
}