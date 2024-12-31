#pragma once

//std
#include <cstdint>

namespace math
{
	class matrix;
}

namespace math
{
	class span
	{
	public:
		//constructors
		span(matrix&, uint32_t, uint32_t, uint32_t = 3, uint32_t = 3);

		//destructor
		~span(void);

		//operators
		const span& operator=(double) const;
		const span& operator=(const span&) const;
		const span& operator=(const matrix&) const;

		const span& operator+=(double) const;
		const span& operator+=(const span&) const;
		const span& operator+=(const matrix&) const;

		const span& operator-=(double) const;
		const span& operator-=(const span&) const;
		const span& operator-=(const matrix&) const;

		double& operator()(uint32_t, uint32_t) const;

		//friends
		friend matrix operator*(const span&, const matrix&);
		friend matrix operator*(const matrix&, const span&);

	protected:
		//data
		matrix& m_k;
		uint32_t m_row;
		uint32_t m_col;
		uint32_t m_rows;
		uint32_t m_cols;
	};
}