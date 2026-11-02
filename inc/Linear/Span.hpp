#pragma once

//std
#include <cstdint>

namespace math
{
	class Mat3;
	class Matrix;
}

namespace math
{
	class Span
	{
	public:
		//constructors
		Span(Matrix&, uint32_t, uint32_t, uint32_t, uint32_t);

		//destructor
		~Span(void) = default;

		//utility
		void eye(void) const;
		void zeros(void) const;

		//operators
		const Span& operator=(double) const;
		const Span& operator=(const Span&) const;
		const Span& operator=(const Mat3&) const;
		const Span& operator=(const Matrix&) const;

		const Span& operator+=(double) const;
		const Span& operator+=(const Span&) const;
		const Span& operator+=(const Mat3&) const;
		const Span& operator+=(const Matrix&) const;

		const Span& operator-=(double) const;
		const Span& operator-=(const Span&) const;
		const Span& operator-=(const Mat3&) const;
		const Span& operator-=(const Matrix&) const;

		double& operator()(uint32_t, uint32_t) const;

		//friends
		friend Matrix operator*(const Span&, const Matrix&);
		friend Matrix operator*(const Matrix&, const Span&);

	protected:
		//data
		Matrix& m_K;
		uint32_t m_row;
		uint32_t m_col;
		uint32_t m_rows;
		uint32_t m_cols;
	};
}