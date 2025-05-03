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
	sparse::sparse(uint32_t rows, uint32_t cols) : 
		m_own(true), m_data_ptr(nullptr), m_data_ref(nullptr), m_rows(rows), m_cols(cols),
		m_rows_map_ptr(nullptr), m_cols_map_ptr(nullptr), m_rows_map_ref(nullptr), m_cols_map_ref(nullptr)
	{
		return;
	}
	sparse::sparse(double* data_ptr, int32_t* rows_map, int32_t* cols_map, uint32_t rows, uint32_t cols) : 
		m_own(false), m_data_ptr(data_ptr), m_data_ref(data_ptr), m_rows(rows), m_cols(cols), 
		m_rows_map_ptr(rows_map), m_cols_map_ptr(cols_map), m_rows_map_ref(rows_map), m_cols_map_ref(cols_map)
	{
		return;
	}
	sparse::sparse(double* data_ptr, const int32_t* rows_map, const int32_t* cols_map, uint32_t rows, uint32_t cols) : 
		m_own(false), m_data_ptr(data_ptr), m_data_ref(data_ptr), m_rows(rows), m_cols(cols), 
		m_rows_map_ptr(nullptr), m_cols_map_ptr(nullptr), m_rows_map_ref(rows_map), m_cols_map_ref(cols_map)
	{
		return;
	}
	sparse::sparse(const double* data_ref, const int32_t* rows_map, const int32_t* cols_map, uint32_t rows, uint32_t cols) : 
		m_own(false), m_data_ptr(nullptr), m_data_ref(data_ref), m_rows(rows), m_cols(cols), 
		m_rows_map_ptr(nullptr), m_cols_map_ptr(nullptr), m_rows_map_ref(rows_map), m_cols_map_ref(cols_map)
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
	void sparse::pattern(int32_t* rows_map, int32_t* cols_map)
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
	void sparse::pattern(int32_t*& rows_map, int32_t*& cols_map) const
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

	bool sparse::solve(vector& x, const vector& f) const
	{
		//check
		const uint32_t n = m_rows;
		int32_t test = !UMFPACK_OK;
		if(m_cols != x.rows() || m_rows != f.rows())
		{
			fprintf(stderr, "Error: Sparse solve called with incompatible vectors!\n");
			exit(EXIT_FAILURE);
		}
		//solve
		double* xd = x.data();
		const double* fd = f.data();
		const double* Kd = m_data_ref;
		const int32_t* r = m_rows_map_ref;
		const int32_t* c = m_cols_map_ref;
		if(umfpack_di_symbolic(n, n, c, r, Kd, &m_symbolic, nullptr, nullptr) == UMFPACK_OK)
		{
			if(umfpack_di_numeric(c, r, Kd, m_symbolic, &m_numeric, nullptr, nullptr) == UMFPACK_OK)
			{
				test = umfpack_di_solve(UMFPACK_A, c, r, Kd, xd, fd, m_numeric, nullptr, nullptr);
			}
		}
		//free memory
		umfpack_di_free_numeric(&m_numeric);
		umfpack_di_free_symbolic(&m_symbolic);
		//return
		return test == UMFPACK_OK;
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
	void sparse::span(sparse& matrix, uint32_t i, uint32_t j) const
	{
		//data
		span_check(matrix, i, j, true);
		const uint32_t nc = matrix.m_cols;
		const uint32_t nr = matrix.m_rows;
		//setup
		span_count(matrix, i, j, true);
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
	void sparse::span_data(sparse& matrix, uint32_t i, uint32_t j) const
	{
		//data
		span_check(matrix, i, j, false);
		const uint32_t nr = matrix.m_rows;
		const uint32_t nc = matrix.m_cols;
		//setup
		for(uint32_t a = 0; a < nc; a++)
		{
			uint32_t counter = 0;
			for(int32_t b = m_cols_map_ref[j + a]; b < m_cols_map_ref[j + a + 1]; b++)
			{
				if(m_rows_map_ref[b] >= int32_t(i) && m_rows_map_ref[b] < int32_t(i + nr))
				{
					matrix.m_data_ptr[matrix.m_cols_map_ref[a] + counter] = m_data_ref[b];
					counter++;
				}
			}
		}
	}
	//sparse
	void sparse::span_pattern(sparse& matrix, uint32_t i, uint32_t j) const
	{
		//data
		span_check(matrix, i, j, true);
		const uint32_t nr = matrix.m_rows;
		const uint32_t nc = matrix.m_cols;
		//setup
		span_count(matrix, i, j, false);
		for(uint32_t a = 0; a < nc; a++)
		{
			uint32_t counter = 0;
			for(int32_t b = m_cols_map_ref[j + a]; b < m_cols_map_ref[j + a + 1]; b++)
			{
				if(m_rows_map_ref[b] >= int32_t(i) && m_rows_map_ref[b] < int32_t(i + nr))
				{
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
			delete[] m_rows_map_ptr;
			delete[] m_cols_map_ptr;
			m_data_ref = m_data_ptr = nullptr;
			m_rows_map_ref = m_rows_map_ptr = nullptr;
			m_cols_map_ref = m_cols_map_ptr = nullptr;
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

	//span
	void sparse::span_count(sparse& matrix, uint32_t i, uint32_t j, bool data) const
	{
		//data
		matrix.cleanup();
		const uint32_t nc = matrix.m_cols;
		const uint32_t nr = matrix.m_rows;
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
		if(data) matrix.m_data_ref = matrix.m_data_ptr = new double[matrix.m_cols_map_ptr[nc]];
		matrix.m_rows_map_ref = matrix.m_rows_map_ptr = new int32_t[matrix.m_cols_map_ptr[nc]];
	}
	void sparse::span_check(sparse& matrix, uint32_t i, uint32_t j, bool owner_check) const
	{
		if(owner_check && !matrix.m_own)
		{
			printf("Error: Sparse matrix span called on a matrix that don't owns data!\n");
			exit(EXIT_FAILURE);
		}
		if(i + matrix.m_rows > m_rows || j + matrix.m_cols > m_cols)
		{
			printf("Error: Sparse matrix span called with uncompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
	}
}