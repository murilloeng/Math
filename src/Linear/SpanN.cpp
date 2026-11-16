//std
#include <stdexcept>

//Math
#include "Math/inc/Linear/SpanN.hpp"
#include "Math/inc/Linear/MatrixN.hpp"
#include "Math/inc/Linear/VectorN.hpp"

namespace math
{
	//constructors
	SpanN::SpanN(MatrixN& K, uint32_t row, uint32_t col, uint32_t rows, uint32_t cols) : m_K{K}, m_row{row}, m_col{col}, m_rows{rows}, m_cols{cols}
	{
		return;
	}

	//linear
	void SpanN::eye(void) const
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				(*this)(i, j) = i == j;
			}
		}
	}
	void SpanN::zeros(void) const
	{
		*this = 0;
	}

	//opertors
	const SpanN& SpanN::operator=(double v) const
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
	const SpanN& SpanN::operator=(const SpanN& s) const
	{
		if(m_rows != s.m_rows || m_cols != s.m_cols)
		{
			throw std::runtime_error("SpanN assign has incompatible dimensions!");
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
	const SpanN& SpanN::operator=(const MatrixN& M) const
	{
		if(m_rows != M.m_rows || m_cols != M.m_cols)
		{
			throw std::runtime_error("SpanN assign has incompatible dimensions!");
		}
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				(*this)(i, j) = M(i, j);
			}
		}
		return *this;
	}
	const SpanN& SpanN::operator=(const VectorN& v) const
	{
		if(m_rows != v.m_rows || m_cols != 1)
		{
			throw std::runtime_error("SpanN assign has incompatible dimensions!");
		}
		for(uint32_t i = 0; i < m_rows; i++)
		{
			(*this)(i, 0) = v[i];
		}
		return *this;
	}

	const SpanN& SpanN::operator+=(double v) const
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
	const SpanN& SpanN::operator+=(const SpanN& s) const
	{
		if(m_rows != s.m_rows || m_cols != s.m_cols)
		{
			throw std::runtime_error("SpanN increment has incompatible dimensions!");
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
	const SpanN& SpanN::operator+=(const MatrixN& k) const
	{
		if(m_rows != k.m_rows || m_cols != k.m_cols)
		{
			throw std::runtime_error("SpanN increment has incompatible dimensions!");
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

	const SpanN& SpanN::operator-=(double v) const
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
	const SpanN& SpanN::operator-=(const SpanN& s) const
	{
		if(s.m_rows != m_rows || s.m_cols != m_cols)
		{
			throw std::runtime_error("SpanN increment has incompatible dimensions!");
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
	const SpanN& SpanN::operator-=(const MatrixN& k) const
	{
		if(m_rows != k.m_rows || m_cols != k.m_cols)
		{
			throw std::runtime_error("SpanN increment has incompatible dimensions!");
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

	double& SpanN::operator()(uint32_t i, uint32_t j) const
	{
		return m_K(i + m_row, j + m_col);
	}

	//friends
	MatrixN operator*(const SpanN& s, const MatrixN& k)
	{
		if(s.m_cols != k.m_rows)
		{
			throw std::runtime_error("SpanN-MatrixN product has incompatible dimensions!");
		}
		MatrixN r(s.m_rows, k.m_cols);
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
	MatrixN operator*(const MatrixN& k, const SpanN& s)
	{
		if(k.m_cols != s.m_rows)
		{
			throw std::runtime_error("SpanN-MatrixN product has incompatible dimensions!");
		}
		MatrixN r(k.m_rows, s.m_cols);
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