#pragma once

//std
#include <cstdint>

namespace math
{
	class MatrixX;
	class VectorX;
}

namespace math
{
	class SpanX
	{
	public:
		//constructors
		SpanX(MatrixX&, uint32_t, uint32_t, uint32_t, uint32_t);

		//destructor
		~SpanX(void) = default;

		//utility
		void eye(void) const;
		void zeros(void) const;

		//operators
		const SpanX& operator=(double) const;
		const SpanX& operator=(const SpanX&) const;
		const SpanX& operator=(const MatrixX&) const;
		const SpanX& operator=(const VectorX&) const;

		const SpanX& operator+=(double) const;
		const SpanX& operator+=(const SpanX&) const;
		const SpanX& operator+=(const MatrixX&) const;
		const SpanX& operator+=(const VectorX&) const;

		const SpanX& operator-=(double) const;
		const SpanX& operator-=(const SpanX&) const;
		const SpanX& operator-=(const MatrixX&) const;
		const SpanX& operator-=(const VectorX&) const;

		double& operator()(uint32_t, uint32_t) const;

		//friends
		friend MatrixX operator*(const SpanX&, const MatrixX&);
		friend MatrixX operator*(const MatrixX&, const SpanX&);

	protected:
		//data
		MatrixX& m_K;
		uint32_t m_row;
		uint32_t m_col;
		uint32_t m_rows;
		uint32_t m_cols;
	};
}