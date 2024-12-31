//std
#include <cstdio>
#include <cstdlib>

//math
#include "Galileo/mat/inc/linear/span.hpp"
#include "Galileo/mat/inc/linear/mat3.hpp"
#include "Galileo/mat/inc/linear/matrix.hpp"

namespace math
{
	//constructors
	span::span(matrix& k, uint32_t row, uint32_t col, uint32_t rows, uint32_t cols) : 
		m_k(k), m_row(row), m_col(col), m_rows(rows), m_cols(cols)
	{
		return;
	}

	//destructor
	span::~span(void)
	{
		return;
	}

	//opertors
	const span& span::operator=(double v) const
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				(*this)(i, j) = v;
			}
		}
		return *this;
	}
	const span& span::operator=(const span& s) const
	{
		if(m_rows != s.m_rows || m_cols != s.m_cols)
		{
			fprintf(stderr, "Error: Span assign with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				(*this)(i, j) = s(i, j);
			}
		}
		return *this;
	}
	const span& span::operator=(const matrix& k) const
	{
		if(m_rows != k.m_rows || m_cols != k.m_cols)
		{
			fprintf(stderr, "Error: Span assign with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				(*this)(i, j) = k(i, j);
			}
		}
		return *this;
	}

	const span& span::operator+=(double v) const
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				(*this)(i, j) += v;
			}
		}
		return *this;
	}
	const span& span::operator+=(const span& s) const
	{
		if(m_rows != s.m_rows || m_cols != s.m_cols)
		{
			fprintf(stderr, "Error: Span increment with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				(*this)(i, j) += s(i, j);
			}
		}
		return *this;
	}
	const span& span::operator+=(const matrix& k) const
	{
		if(m_rows != k.m_rows || m_cols != k.m_cols)
		{
			fprintf(stderr, "Error: Span increment with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				(*this)(i, j) += k(i, j);
			}
		}
		return *this;
	}

	const span& span::operator-=(double v) const
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				(*this)(i, j) -= v;
			}
		}
		return *this;
	}
	const span& span::operator-=(const span& s) const
	{
		if(s.m_rows != m_rows || s.m_cols != m_cols)
		{
			fprintf(stderr, "Error: Span increment with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				(*this)(i, j) -= s(i, j);
			}
		}
		return *this;
	}
	const span& span::operator-=(const matrix& k) const
	{
		if(m_rows != k.m_rows || m_cols != k.m_cols)
		{
			fprintf(stderr, "Error: Span increment with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				(*this)(i, j) -= k(i, j);
			}
		}
		return *this;
	}

	double& span::operator()(uint32_t i, uint32_t j) const
	{
		return m_k(i + m_row, j + m_col);
	}

	//friends
	matrix operator*(const span& s, const matrix& k)
	{
		if(s.m_cols != k.m_rows)
		{
			fprintf(stderr, "Error: Span-Matrix product with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		matrix r(s.m_rows, k.m_cols);
		for(uint32_t i = 0; i < s.m_rows; i++)
		{
			for(uint32_t j = 0; j < k.m_cols; j++)
			{
				r(i, j) = 0;
				for(uint32_t p = 0; p < s.m_cols; p++)
				{
					r(i, j) += s(i, p) * k(p, j);
				}
			}
		}
		return r;
	}
	matrix operator*(const matrix& k, const span& s)
	{
		if(k.m_cols != s.m_rows)
		{
			fprintf(stderr, "Error: Matrix-Span product with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		matrix r(k.m_rows, s.m_cols);
		for(uint32_t i = 0; i < k.m_rows; i++)
		{
			for(uint32_t j = 0; j < s.m_cols; j++)
			{
				r(i, j) = 0;
				for(uint32_t p = 0; p < k.m_cols; p++)
				{
					r(i, j) += k(i, p) * s(p, j);
				}
			}
		}
		return r;
	}
}