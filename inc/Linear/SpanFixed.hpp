#pragma once

//std
#include <cstdint>
#include <stdexcept>

//Math
#include "Math/inc/Linear/Matrix.hpp"

namespace math
{
	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c) class SpanFixed
	{
	public:
		//constructors
		SpanFixed(Matrix<r, c>&, uint32_t, uint32_t);
		SpanFixed(const SpanFixed&) = default;

		//destructor
		~SpanFixed(void) = default;

		//utility
		void eye(void) const requires (m == n);
		void zeros(void) const;

		//operators
		const SpanFixed& operator=(double) const;
		const SpanFixed& operator=(const SpanFixed&) const;
		const SpanFixed& operator=(const Matrix<m, n>&) const;
		template<uint32_t a, uint32_t b> const SpanFixed& operator=(const SpanFixed<a, b, m, n>&) const;

		const SpanFixed& operator+=(double) const;
		const SpanFixed& operator+=(const Matrix<m, n>&) const;
		template<uint32_t a, uint32_t b> const SpanFixed& operator+=(const SpanFixed<a, b, m, n>&) const;

		const SpanFixed& operator-=(double) const;
		const SpanFixed& operator-=(const Matrix<m, n>&) const;
		template<uint32_t a, uint32_t b> const SpanFixed& operator-=(const SpanFixed<a, b, m, n>&) const;

		double& operator()(uint32_t, uint32_t) const;

		//products
		template<uint32_t p> Matrix<m, p> operator*(const Matrix<n, p>&) const;

	protected:
		//data
		uint32_t m_row;
		uint32_t m_col;
		Matrix<r, c>& m_K;
	};

	//constructors
	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c)
	SpanFixed<r, c, m, n>::SpanFixed(Matrix<r, c>& K, uint32_t row, uint32_t col) : m_row{row}, m_col{col}, m_K{K}
	{
		if(row > r - m || col > c - n)
		{
			throw std::out_of_range("SpanFixed exceeds matrix bounds!");
		}
	}

	//utility
	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c)
	void SpanFixed<r, c, m, n>::eye(void) const requires (m == n)
	{
		for(uint32_t i = 0; i < m; i++)
		{
			for(uint32_t j = 0; j < n; j++)
			{
				(*this)(i, j) = i == j;
			}
		}
	}
	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c)
	void SpanFixed<r, c, m, n>::zeros(void) const
	{
		*this = 0;
	}

	//operators
	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c)
	const SpanFixed<r, c, m, n>& SpanFixed<r, c, m, n>::operator=(double s) const
	{
		for(uint32_t i = 0; i < m; i++)
		{
			for(uint32_t j = 0; j < n; j++)
			{
				(*this)(i, j) = s;
			}
		}
		return *this;
	}
	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c)
	const SpanFixed<r, c, m, n>& SpanFixed<r, c, m, n>::operator=(const SpanFixed& S) const
	{
		return operator=<r, c>(S);
	}
	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c)
	const SpanFixed<r, c, m, n>& SpanFixed<r, c, m, n>::operator=(const Matrix<m, n>& M) const
	{
		for(uint32_t i = 0; i < m; i++)
		{
			for(uint32_t j = 0; j < n; j++)
			{
				(*this)(i, j) = M(i, j);
			}
		}
		return *this;
	}
	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c)
	template<uint32_t a, uint32_t b> const SpanFixed<r, c, m, n>& SpanFixed<r, c, m, n>::operator=(const SpanFixed<a, b, m, n>& V) const
	{
		for(uint32_t i = 0; i < m; i++)
		{
			for(uint32_t j = 0; j < n; j++)
			{
				(*this)(i, j) = V(i, j);
			}
		}
		return *this;
	}

	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c)
	const SpanFixed<r, c, m, n>& SpanFixed<r, c, m, n>::operator+=(double s) const
	{
		for(uint32_t i = 0; i < m; i++)
		{
			for(uint32_t j = 0; j < n; j++)
			{
				(*this)(i, j) += s;
			}
		}
		return *this;
	}
	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c)
	const SpanFixed<r, c, m, n>& SpanFixed<r, c, m, n>::operator+=(const Matrix<m, n>& M) const
	{
		for(uint32_t i = 0; i < m; i++)
		{
			for(uint32_t j = 0; j < n; j++)
			{
				(*this)(i, j) += M(i, j);
			}
		}
		return *this;
	}
	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c)
	template<uint32_t a, uint32_t b> const SpanFixed<r, c, m, n>& SpanFixed<r, c, m, n>::operator+=(const SpanFixed<a, b, m, n>& V) const
	{
		for(uint32_t i = 0; i < m; i++)
		{
			for(uint32_t j = 0; j < n; j++)
			{
				(*this)(i, j) += V(i, j);
			}
		}
		return *this;
	}

	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c)
	const SpanFixed<r, c, m, n>& SpanFixed<r, c, m, n>::operator-=(double s) const
	{
		for(uint32_t i = 0; i < m; i++)
		{
			for(uint32_t j = 0; j < n; j++)
			{
				(*this)(i, j) -= s;
			}
		}
		return *this;
	}
	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c)
	const SpanFixed<r, c, m, n>& SpanFixed<r, c, m, n>::operator-=(const Matrix<m, n>& M) const
	{
		for(uint32_t i = 0; i < m; i++)
		{
			for(uint32_t j = 0; j < n; j++)
			{
				(*this)(i, j) -= M(i, j);
			}
		}
		return *this;
	}
	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c)
	template<uint32_t a, uint32_t b> const SpanFixed<r, c, m, n>& SpanFixed<r, c, m, n>::operator-=(const SpanFixed<a, b, m, n>& V) const
	{
		for(uint32_t i = 0; i < m; i++)
		{
			for(uint32_t j = 0; j < n; j++)
			{
				(*this)(i, j) -= V(i, j);
			}
		}
		return *this;
	}

	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c)
	double& SpanFixed<r, c, m, n>::operator()(uint32_t i, uint32_t j) const
	{
		return m_K(m_row + i, m_col + j);
	}

	//products
	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c)
	template<uint32_t p> Matrix<m, p> SpanFixed<r, c, m, n>::operator*(const Matrix<n, p>& M) const
	{
		Matrix<m, p> A;
		for(uint32_t i = 0; i < m; i++)
		{
			for(uint32_t j = 0; j < p; j++)
			{
				A(i, j) = 0;
				for(uint32_t k = 0; k < n; k++)
				{
					A(i, j) += (*this)(i, k) * M(k, j);
				}
			}
		}
		return A;
	}
	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n, uint32_t p> requires (m <= r && n <= c)
	Matrix<p, n> operator*(const Matrix<p, m>& M, const SpanFixed<r, c, m, n>& S)
	{
		Matrix<p, n> A;
		for(uint32_t i = 0; i < p; i++)
		{
			for(uint32_t j = 0; j < n; j++)
			{
				A(i, j) = 0;
				for(uint32_t k = 0; k < m; k++)
				{
					A(i, j) += M(i, k) * S(k, j);
				}
			}
		}
		return A;
	}
}
