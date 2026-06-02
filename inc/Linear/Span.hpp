#pragma once

//std
#include <cstdint>

namespace math
{
	class Matrix;
}

namespace math
{
	class Span
	{
	public:
		//constructors
		Span(Matrix&, uint32_t, uint32_t, uint32_t = 3, uint32_t = 3);

		//destructor
		~Span(void);

		//operators
		const Span& operator=(double) const;
		const Span& operator=(const Span&) const;
		const Span& operator=(const Matrix&) const;

		const Span& operator+=(double) const;
		const Span& operator+=(const Span&) const;
		const Span& operator+=(const Matrix&) const;

		const Span& operator-=(double) const;
		const Span& operator-=(const Span&) const;
		const Span& operator-=(const Matrix&) const;

		double& operator()(uint32_t, uint32_t) const;

		//friends
		friend Matrix operator*(const Span&, const Matrix&);
		friend Matrix operator*(const Matrix&, const Span&);

	protected:
		//data
		Matrix& m_k;
		uint32_t m_row;
		uint32_t m_col;
		uint32_t m_rows;
		uint32_t m_cols;
	};
}