#pragma once

//std
#include <cstdint>

namespace math
{
	class vector;
	class matrix;
}

namespace math
{
	class sparse
	{
	public:
		//constructors
		sparse(double*, const int32_t*, const int32_t*, uint32_t, uint32_t);
		sparse(const double*, const int32_t*, const int32_t*, uint32_t, uint32_t);

		//destructor
		~sparse(void);

		//data
		double* data(void);
		const double* data(void) const;

		//linear
		double norm(void) const;
		double trace(void) const;

		//print
		void print(const char* = "", bool = true) const;

		//operators
		vector operator*(const vector&) const;
		double& operator()(uint32_t, uint32_t);
		const double& operator()(uint32_t, uint32_t) const;

		//convert
		matrix convert(void) const;

	private:
		//data
		uint32_t m_rows;
		uint32_t m_cols;

		double* m_ptr;
		const double* m_ref;

		const int32_t* m_row_map;
		const int32_t* m_col_map;
	};
}