#pragma once

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
		span(matrix&, unsigned, unsigned, unsigned = 3, unsigned = 3);

		//destructor
		virtual ~span(void);

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

		double& operator()(unsigned, unsigned) const;

		//friends
		friend matrix operator*(const span&, const matrix&);
		friend matrix operator*(const matrix&, const span&);

	protected:
		//data
		matrix& m_k;
		unsigned m_row;
		unsigned m_col;
		unsigned m_rows;
		unsigned m_cols;
	};
}