//std
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

//Math
#include "Math/inc/Linear/Span.hpp"
#include "Math/inc/Linear/Mat3.hpp"
#include "Math/inc/Linear/Matrix.hpp"

namespace math
{
	//constructors
	Span::Span(Matrix& k, uint32_t row, uint32_t col, uint32_t rows, uint32_t cols) : 
		m_k(k), m_row(row), m_col(col), m_rows(rows), m_cols(cols)
	{
		return;
	}

	//destructor
	Span::~Span(void)
	{
		return;
	}

	//opertors
	const Span& Span::operator=(double v) const
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
	const Span& Span::operator=(const Span& s) const
	{
		if(m_rows != s.m_rows || m_cols != s.m_cols)
		{
			throw std::runtime_error("Span assign has incompatible dimensions!");
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
	const Span& Span::operator=(const Matrix& k) const
	{
		if(m_rows != k.m_rows || m_cols != k.m_cols)
		{
			throw std::runtime_error("Span assign has incompatible dimensions!");
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

	const Span& Span::operator+=(double v) const
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
	const Span& Span::operator+=(const Span& s) const
	{
		if(m_rows != s.m_rows || m_cols != s.m_cols)
		{
			throw std::runtime_error("Span increment has incompatible dimensions!");
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
	const Span& Span::operator+=(const Matrix& k) const
	{
		if(m_rows != k.m_rows || m_cols != k.m_cols)
		{
			throw std::runtime_error("Span increment has incompatible dimensions!");
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

	const Span& Span::operator-=(double v) const
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
	const Span& Span::operator-=(const Span& s) const
	{
		if(s.m_rows != m_rows || s.m_cols != m_cols)
		{
			throw std::runtime_error("Span increment has incompatible dimensions!");
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
	const Span& Span::operator-=(const Matrix& k) const
	{
		if(m_rows != k.m_rows || m_cols != k.m_cols)
		{
			throw std::runtime_error("Span increment has incompatible dimensions!");
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

	double& Span::operator()(uint32_t i, uint32_t j) const
	{
		return m_k(i + m_row, j + m_col);
	}

	//friends
	Matrix operator*(const Span& s, const Matrix& k)
	{
		if(s.m_cols != k.m_rows)
		{
			throw std::runtime_error("Span-Matrix product has incompatible dimensions!");
		}
		Matrix r(s.m_rows, k.m_cols);
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
	Matrix operator*(const Matrix& k, const Span& s)
	{
		if(k.m_cols != s.m_rows)
		{
			throw std::runtime_error("Span-Matrix product has incompatible dimensions!");
		}
		Matrix r(k.m_rows, s.m_cols);
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