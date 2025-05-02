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
	sparse::sparse(uint32_t rows, uint32_t cols) : m_own(true), 
		m_data_ptr(nullptr), m_data_ref(nullptr), m_rows(rows), m_cols(cols),
		m_rows_map_ptr(nullptr), m_cols_map_ptr(nullptr), m_rows_map_ref(nullptr), m_cols_map_ref(nullptr)
	{
		return;
	}
	sparse::sparse(double* ptr, int32_t* row_map, int32_t* col_map, uint32_t rows, uint32_t cols) : 
		m_own(false), m_data_ptr(ptr), m_data_ref(ptr), m_rows(rows), m_cols(cols), 
		m_rows_map_ptr(row_map), m_cols_map_ptr(col_map), m_rows_map_ref(row_map), m_cols_map_ref(col_map)
	{
		return;
	}
	sparse::sparse(const double* ref, const int32_t* row_map, const int32_t* col_map, uint32_t rows, uint32_t cols) : 
		m_own(false), m_data_ptr(nullptr), m_data_ref(ref), m_rows(rows), m_cols(cols), 
		m_rows_map_ptr(nullptr), m_cols_map_ptr(nullptr), m_rows_map_ref(row_map), m_cols_map_ref(col_map)
	{
		return;
	}

	//constructors
	sparse::~sparse(void)
	{
		cleanup();
	}

	//data
	double* sparse::data(void)
	{
		return m_data_ptr;
	}
	const double* sparse::data(void) const
	{
		return m_data_ref;
	}

	uint32_t sparse::rows(void) const
	{
		return m_rows;
	}
	uint32_t sparse::cols(void) const
	{
		return m_cols;
	}

	int32_t* sparse::rows_map(void)
	{
		return m_rows_map_ptr;
	}
	int32_t* sparse::cols_map(void)
	{
		return m_cols_map_ptr;
	}
	const int32_t* sparse::rows_map(void) const
	{
		return m_rows_map_ref;
	}
	const int32_t* sparse::cols_map(void) const
	{
		return m_cols_map_ref;
	}

	//pattern
	void sparse::pattern(int32_t* cols_map, int32_t* rows_map)
	{
		//check
		if(!m_own)
		{
			printf("Error: Sparse matrix pattern change called on an object that don owns its data!\n");
			exit(EXIT_FAILURE);
		}
		//setup
		cleanup();
		m_data_ref = m_data_ptr = new double[cols_map[m_cols]];
		m_cols_map_ref = m_cols_map_ptr = new int32_t[m_cols + 1];
		m_rows_map_ref = m_rows_map_ptr = new int32_t[cols_map[m_cols]];
		//pattern
		memcpy(m_cols_map_ptr, cols_map, (m_cols + 1) * sizeof(int32_t));
		memcpy(m_rows_map_ptr, rows_map, cols_map[m_cols] * sizeof(int32_t));
	}
	void sparse::pattern(int32_t*& cols_map, int32_t*& rows_map) const
	{
		delete[] cols_map;
		delete[] rows_map;
		cols_map = new int32_t[m_cols + 1];
		rows_map = new int32_t[m_cols_map_ref[m_cols]];
		memcpy(cols_map, m_cols_map_ref, (m_cols + 1) * sizeof(int32_t));
		memcpy(rows_map, m_rows_map_ref, m_cols_map_ref[m_cols] * sizeof(int32_t));
	}

	//linear
	double sparse::norm(void) const
	{
		double s = 0;
		for(uint32_t i = 0; i < m_cols; i++)
		{
			for(int32_t j = m_cols_map_ref[i]; j < m_cols_map_ref[i + 1]; j++)
			{
				s += m_data_ref[j] * m_data_ref[j];
			}
		}
		return sqrt(s);
	}
	double sparse::trace(void) const
	{
		double s = 0;
		for(uint32_t i = 0; i < m_cols; i++)
		{
			for(int32_t j = m_cols_map_ref[i]; j < m_cols_map_ref[i + 1]; j++)
			{
				if(i == (uint32_t) m_rows_map_ref[j])
				{
					s += m_data_ref[j];
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
		dense ? print_dense() : print_sparse();
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
			for(int32_t j = m_cols_map_ref[i]; j < m_cols_map_ref[i + 1]; j++)
			{
				r[m_rows_map_ref[j]] += m_data_ref[j] * v[i];
			}
		}
		return r;
	}
	double& sparse::operator()(uint32_t i, uint32_t j)
	{
		for(int32_t k = m_cols_map_ref[j]; k < m_cols_map_ref[j + 1]; k++)
		{
			if(int32_t(i) == m_rows_map_ref[k])
			{
				return m_data_ptr[k];
			}
		}
		fprintf(stderr, "Error: Sparse matrix operator() called with out of range index!\n");
		exit(EXIT_FAILURE);
	}
	const double& sparse::operator()(uint32_t i, uint32_t j) const
	{
		for(int32_t p = m_cols_map_ref[j]; p < m_cols_map_ref[j + 1]; p++)
		{
			if(int32_t(i) == m_rows_map_ref[p])
			{
				return m_data_ref[p];
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
			for(int32_t j = m_cols_map_ref[i]; j < m_cols_map_ref[i + 1]; j++)
			{
				M(m_rows_map_ref[j], i) = m_data_ref[j];
			}
		}
		return M;
	}

	//sparse
	void sparse::span(sparse& matrix, uint32_t i, uint32_t j, uint32_t nr, uint32_t nc) const
	{
		//check
		if(!matrix.m_own)
		{
			printf("Error: Sparse matrix span called on a matrix that don't owns data!\n");
			exit(EXIT_FAILURE);
		}
		if(i + nr > m_rows || j + nc > m_cols)
		{
			printf("Error: Sparse matrix span called with uncompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		//setup
		matrix.cleanup();
		matrix.m_cols_map_ref = matrix.m_cols_map_ptr = new int32_t[nc + 1];
		//count
		matrix.m_cols_map_ptr[0] = 0;
		for(uint32_t a = 0; a < nc; a++)
		{
			matrix.m_cols_map_ptr[a + 1] = matrix.m_cols_map_ptr[a];
			for(int32_t b = m_cols_map_ref[j + a]; b < m_cols_map_ref[j + a + 1]; b++)
			{
				if(m_rows_map_ref[b] >= int32_t(i) && m_rows_map_ref[b] < int32_t(i + nr))
				{
					matrix.m_cols_map_ptr[a + 1]++;
				}
			}
		}
		matrix.m_data_ref = matrix.m_data_ptr = new double[matrix.m_cols_map_ptr[nc]];
		matrix.m_rows_map_ref = matrix.m_rows_map_ptr = new int32_t[matrix.m_cols_map_ptr[nc]];
		//setup
		for(uint32_t a = 0; a < nc; a++)
		{
			uint32_t counter = 0;
			for(int32_t b = m_cols_map_ref[j + a]; b < m_cols_map_ref[j + a + 1]; b++)
			{
				if(m_rows_map_ref[b] >= int32_t(i) && m_rows_map_ref[b] < int32_t(i + nr))
				{
					matrix.m_data_ptr[matrix.m_cols_map_ptr[a] + counter] = m_data_ref[b];
					matrix.m_rows_map_ptr[matrix.m_cols_map_ptr[a] + counter] = m_rows_map_ref[b] - i;
					counter++;
				}
			}
		}
	}
	//sparse
	void sparse::span_data(sparse& matrix, uint32_t i, uint32_t j, uint32_t nr, uint32_t nc) const
	{
		//check
		if(i + nr > m_rows || j + nc > m_cols)
		{
			printf("Error: Sparse matrix span called with uncompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		//setup
		for(uint32_t a = 0; a < nc; a++)
		{
			uint32_t counter = 0;
			for(int32_t b = m_cols_map_ref[j + a]; b < m_cols_map_ref[j + a + 1]; b++)
			{
				if(m_rows_map_ref[b] >= int32_t(i) && m_rows_map_ref[b] < int32_t(i + nr))
				{
					matrix.m_data_ptr[matrix.m_cols_map_ptr[a] + counter] = m_data_ref[b];
					counter++;
				}
			}
		}
	}
	//sparse
	void sparse::span_pattern(sparse& matrix, uint32_t i, uint32_t j, uint32_t nr, uint32_t nc) const
	{
		//check
		if(!matrix.m_own)
		{
			printf("Error: Sparse matrix span called on a matrix that don't owns data!\n");
			exit(EXIT_FAILURE);
		}
		if(i + nr > m_rows || j + nc > m_cols)
		{
			printf("Error: Sparse matrix span called with uncompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		//setup
		delete[] matrix.m_data_ptr;
		delete[] matrix.m_cols_map_ptr;
		delete[] matrix.m_rows_map_ptr;
		matrix.m_cols_map_ref = matrix.m_cols_map_ptr = new int32_t[nc + 1];
		//count
		matrix.m_cols_map_ptr[0] = 0;
		for(uint32_t a = 0; a < nc; a++)
		{
			matrix.m_cols_map_ptr[a + 1] = matrix.m_cols_map_ptr[a];
			for(int32_t b = m_cols_map_ref[j + a]; b < m_cols_map_ref[j + a + 1]; b++)
			{
				if(m_rows_map_ref[b] >= int32_t(i) && m_rows_map_ref[b] < int32_t(i + nr))
				{
					matrix.m_cols_map_ptr[a + 1]++;
				}
			}
		}
		matrix.m_data_ref = matrix.m_data_ptr = new double[matrix.m_cols_map_ptr[nc]];
		matrix.m_rows_map_ref = matrix.m_rows_map_ptr = new int32_t[matrix.m_cols_map_ptr[nc]];
		//setup
		for(uint32_t a = 0; a < nc; a++)
		{
			uint32_t counter = 0;
			for(int32_t b = m_cols_map_ref[j + a]; b < m_cols_map_ref[j + a + 1]; b++)
			{
				if(m_rows_map_ref[b] >= int32_t(i) && m_rows_map_ref[b] < int32_t(i + nr))
				{
					matrix.m_data_ptr[matrix.m_cols_map_ptr[a] + counter] = m_data_ref[b];
					matrix.m_rows_map_ptr[matrix.m_cols_map_ptr[a] + counter] = m_rows_map_ref[b] - i;
					counter++;
				}
			}
		}
	}

	//data
	void sparse::cleanup(void)
	{
		if(m_own)
		{
			delete[] m_data_ptr;
			delete[] m_data_ref;
			delete[] m_rows_map_ptr;
			delete[] m_rows_map_ptr;
			delete[] m_cols_map_ptr;
			delete[] m_cols_map_ref;
		}
	}

	//print
	void sparse::print_dense(void) const
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
	void sparse::print_sparse(void) const
	{
		for(uint32_t i = 0; i < m_cols; i++)
		{
			for(int32_t j = m_cols_map_ref[i]; j < m_cols_map_ref[i + 1]; j++)
			{
				printf("(%04d, %04d): %+.2e\n", m_rows_map_ref[j], i, m_data_ref[j]);
			}
		}
	}
}