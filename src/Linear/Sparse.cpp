//std
#include <cmath>
#include <cstdio>
#include <cstring>
#include <malloc.h>
#include <stdexcept>

//umfpack
#include <suitesparse/umfpack.h>

//Math
#include "Math/inc/Linear/Sparse.hpp"
#include "Math/inc/Linear/Vector.hpp"

namespace math
{
	//constructors
	Sparse::Sparse(uint32_t rows, uint32_t cols) : 
		m_own(true), m_data_ptr(nullptr), m_data_ref(nullptr), m_rows(rows), m_cols(cols),
		m_rows_map_ptr(nullptr), m_cols_map_ptr(nullptr), m_rows_map_ref(nullptr), m_cols_map_ref(nullptr)
	{
		return;
	}
	Sparse::Sparse(double* data_ptr, int32_t* rows_map, int32_t* cols_map, uint32_t rows, uint32_t cols) : 
		m_own(false), m_data_ptr(data_ptr), m_data_ref(data_ptr), m_rows(rows), m_cols(cols), 
		m_rows_map_ptr(rows_map), m_cols_map_ptr(cols_map), m_rows_map_ref(rows_map), m_cols_map_ref(cols_map)
	{
		return;
	}
	Sparse::Sparse(double* data_ptr, const int32_t* rows_map, const int32_t* cols_map, uint32_t rows, uint32_t cols) : 
		m_own(false), m_data_ptr(data_ptr), m_data_ref(data_ptr), m_rows(rows), m_cols(cols), 
		m_rows_map_ptr(nullptr), m_cols_map_ptr(nullptr), m_rows_map_ref(rows_map), m_cols_map_ref(cols_map)
	{
		return;
	}
	Sparse::Sparse(const double* data_ref, const int32_t* rows_map, const int32_t* cols_map, uint32_t rows, uint32_t cols) : 
		m_own(false), m_data_ptr(nullptr), m_data_ref(data_ref), m_rows(rows), m_cols(cols), 
		m_rows_map_ptr(nullptr), m_cols_map_ptr(nullptr), m_rows_map_ref(rows_map), m_cols_map_ref(cols_map)
	{
		return;
	}

	//constructors
	Sparse::~Sparse(void)
	{
		cleanup();
	}

	//data
	double* Sparse::data(void)
	{
		return m_data_ptr;
	}
	const double* Sparse::data(void) const
	{
		return m_data_ref;
	}

	uint32_t Sparse::rows(void) const
	{
		return m_rows;
	}
	uint32_t Sparse::cols(void) const
	{
		return m_cols;
	}

	int32_t* Sparse::rows_map(void)
	{
		return m_rows_map_ptr;
	}
	int32_t* Sparse::cols_map(void)
	{
		return m_cols_map_ptr;
	}
	const int32_t* Sparse::rows_map(void) const
	{
		return m_rows_map_ref;
	}
	const int32_t* Sparse::cols_map(void) const
	{
		return m_cols_map_ref;
	}

	//pattern
	void Sparse::pattern(int32_t* rows_map, int32_t* cols_map)
	{
		//check
		if(!m_own)
		{
			throw std::runtime_error("Sparse pattern change called on Matrix that does not own its memory!");
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
	void Sparse::pattern(int32_t*& rows_map, int32_t*& cols_map) const
	{
		delete[] cols_map;
		delete[] rows_map;
		cols_map = new int32_t[m_cols + 1];
		rows_map = new int32_t[m_cols_map_ref[m_cols]];
		memcpy(cols_map, m_cols_map_ref, (m_cols + 1) * sizeof(int32_t));
		memcpy(rows_map, m_rows_map_ref, m_cols_map_ref[m_cols] * sizeof(int32_t));
	}

	//linear
	double Sparse::norm(void) const
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
	double Sparse::trace(void) const
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

	bool Sparse::solve(Vector& x, const Vector& f) const
	{
		//check
		const uint32_t n = m_rows;
		int32_t test = !UMFPACK_OK;
		if(m_cols != x.rows() || m_rows != f.rows())
		{
			throw std::runtime_error("Sparse solve has incompatible vectors!");
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
	void Sparse::print(const char* header, bool dense) const
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
	Vector Sparse::operator*(const Vector& v) const
	{
		//check
		if(m_cols != v.rows())
		{
			throw std::runtime_error("Sparse Vector product has incompatible dimentions!");
		}
		//vector
		Vector r(v.rows(), mode::zeros);
		for(uint32_t i = 0; i < m_cols; i++)
		{
			for(int32_t j = m_cols_map_ref[i]; j < m_cols_map_ref[i + 1]; j++)
			{
				r[m_rows_map_ref[j]] += m_data_ref[j] * v[i];
			}
		}
		return r;
	}
	double& Sparse::operator()(uint32_t i, uint32_t j)
	{
		for(int32_t k = m_cols_map_ref[j]; k < m_cols_map_ref[j + 1]; k++)
		{
			if(int32_t(i) == m_rows_map_ref[k])
			{
				return m_data_ptr[k];
			}
		}
		throw std::runtime_error("Sparse operator() has index out of range!");
	}
	const double Sparse::operator()(uint32_t i, uint32_t j) const
	{
		for(int32_t p = m_cols_map_ref[j]; p < m_cols_map_ref[j + 1]; p++)
		{
			if(int32_t(i) == m_rows_map_ref[p])
			{
				return m_data_ref[p];
			}
		}
		return 0;
	}

	//convert
	Matrix Sparse::convert(void) const
	{
		Matrix M(m_rows, m_cols, mode::zeros);
		for(uint32_t i = 0; i < m_cols; i++)
		{
			for(int32_t j = m_cols_map_ref[i]; j < m_cols_map_ref[i + 1]; j++)
			{
				M(m_rows_map_ref[j], i) = m_data_ref[j];
			}
		}
		return M;
	}

	//Sparse
	void Sparse::Span(Sparse& Matrix, uint32_t i, uint32_t j) const
	{
		//data
		span_check(Matrix, i, j, true);
		const uint32_t nc = Matrix.m_cols;
		const uint32_t nr = Matrix.m_rows;
		//setup
		span_count(Matrix, i, j, true);
		for(uint32_t a = 0; a < nc; a++)
		{
			uint32_t counter = 0;
			for(int32_t b = m_cols_map_ref[j + a]; b < m_cols_map_ref[j + a + 1]; b++)
			{
				if(m_rows_map_ref[b] >= int32_t(i) && m_rows_map_ref[b] < int32_t(i + nr))
				{
					Matrix.m_data_ptr[Matrix.m_cols_map_ptr[a] + counter] = m_data_ref[b];
					Matrix.m_rows_map_ptr[Matrix.m_cols_map_ptr[a] + counter] = m_rows_map_ref[b] - i;
					counter++;
				}
			}
		}
	}
	//Sparse
	void Sparse::span_data(Sparse& Matrix, uint32_t i, uint32_t j) const
	{
		//data
		span_check(Matrix, i, j, false);
		const uint32_t nr = Matrix.m_rows;
		const uint32_t nc = Matrix.m_cols;
		//setup
		for(uint32_t a = 0; a < nc; a++)
		{
			uint32_t counter = 0;
			for(int32_t b = m_cols_map_ref[j + a]; b < m_cols_map_ref[j + a + 1]; b++)
			{
				if(m_rows_map_ref[b] >= int32_t(i) && m_rows_map_ref[b] < int32_t(i + nr))
				{
					Matrix.m_data_ptr[Matrix.m_cols_map_ref[a] + counter] = m_data_ref[b];
					counter++;
				}
			}
		}
	}
	//Sparse
	void Sparse::span_pattern(Sparse& Matrix, uint32_t i, uint32_t j) const
	{
		//data
		span_check(Matrix, i, j, true);
		const uint32_t nr = Matrix.m_rows;
		const uint32_t nc = Matrix.m_cols;
		//setup
		span_count(Matrix, i, j, false);
		for(uint32_t a = 0; a < nc; a++)
		{
			uint32_t counter = 0;
			for(int32_t b = m_cols_map_ref[j + a]; b < m_cols_map_ref[j + a + 1]; b++)
			{
				if(m_rows_map_ref[b] >= int32_t(i) && m_rows_map_ref[b] < int32_t(i + nr))
				{
					Matrix.m_rows_map_ptr[Matrix.m_cols_map_ptr[a] + counter] = m_rows_map_ref[b] - i;
					counter++;
				}
			}
		}
	}

	//data
	void Sparse::cleanup(void)
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
	void Sparse::print_dense(void) const
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
	void Sparse::print_sparse(void) const
	{
		for(uint32_t i = 0; i < m_cols; i++)
		{
			for(int32_t j = m_cols_map_ref[i]; j < m_cols_map_ref[i + 1]; j++)
			{
				printf("(%04d, %04d): %+.2e\n", m_rows_map_ref[j], i, m_data_ref[j]);
			}
		}
	}

	//Span
	void Sparse::span_count(Sparse& Matrix, uint32_t i, uint32_t j, bool data) const
	{
		//data
		Matrix.cleanup();
		const uint32_t nc = Matrix.m_cols;
		const uint32_t nr = Matrix.m_rows;
		Matrix.m_cols_map_ref = Matrix.m_cols_map_ptr = new int32_t[nc + 1];
		//count
		Matrix.m_cols_map_ptr[0] = 0;
		for(uint32_t a = 0; a < nc; a++)
		{
			Matrix.m_cols_map_ptr[a + 1] = Matrix.m_cols_map_ptr[a];
			for(int32_t b = m_cols_map_ref[j + a]; b < m_cols_map_ref[j + a + 1]; b++)
			{
				if(m_rows_map_ref[b] >= int32_t(i) && m_rows_map_ref[b] < int32_t(i + nr))
				{
					Matrix.m_cols_map_ptr[a + 1]++;
				}
			}
		}
		if(data) Matrix.m_data_ref = Matrix.m_data_ptr = new double[Matrix.m_cols_map_ptr[nc]];
		Matrix.m_rows_map_ref = Matrix.m_rows_map_ptr = new int32_t[Matrix.m_cols_map_ptr[nc]];
	}
	void Sparse::span_check(Sparse& Matrix, uint32_t i, uint32_t j, bool owner_check) const
	{
		if(owner_check && !Matrix.m_own)
		{
			throw std::runtime_error("Sparse Span called on Matrix that does not own its memory!");
		}
		if(i + Matrix.m_rows > m_rows || j + Matrix.m_cols > m_cols)
		{
			throw std::runtime_error("Sparse Span has incompatible dimensions!");
		}
	}
}