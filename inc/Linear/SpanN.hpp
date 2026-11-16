#pragma once

//std
#include <cstdint>

namespace math
{
	class MatrixN;
	class VectorN;
}

namespace math
{
	class SpanN
	{
	public:
		//constructors
		SpanN(MatrixN&, uint32_t, uint32_t, uint32_t, uint32_t);

		//destructor
		~SpanN(void) = default;

		//utility
		void eye(void) const;
		void zeros(void) const;

		//operators
		const SpanN& operator=(double) const;
		const SpanN& operator=(const SpanN&) const;
		const SpanN& operator=(const MatrixN&) const;
		const SpanN& operator=(const VectorN&) const;

		const SpanN& operator+=(double) const;
		const SpanN& operator+=(const SpanN&) const;
		const SpanN& operator+=(const MatrixN&) const;
		const SpanN& operator+=(const VectorN&) const;

		const SpanN& operator-=(double) const;
		const SpanN& operator-=(const SpanN&) const;
		const SpanN& operator-=(const MatrixN&) const;
		const SpanN& operator-=(const VectorN&) const;

		double& operator()(uint32_t, uint32_t) const;

		//friends
		friend MatrixN operator*(const SpanN&, const MatrixN&);
		friend MatrixN operator*(const MatrixN&, const SpanN&);

	protected:
		//data
		MatrixN& m_K;
		uint32_t m_row;
		uint32_t m_col;
		uint32_t m_rows;
		uint32_t m_cols;
	};
}