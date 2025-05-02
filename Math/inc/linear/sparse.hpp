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
		sparse(uint32_t, uint32_t);
		sparse(double*, int32_t*, int32_t*, uint32_t, uint32_t);
		sparse(const double*, const int32_t*, const int32_t*, uint32_t, uint32_t);

		//destructor
		~sparse(void);

		//data
		double* data(void);
		const double* data(void) const;

		uint32_t rows(void) const;
		uint32_t cols(void) const;

		int32_t* rows_map(void);
		int32_t* cols_map(void);
		const int32_t* rows_map(void) const;
		const int32_t* cols_map(void) const;

		//pattern
		void pattern(int32_t*, int32_t*);
		void pattern(int32_t*&, int32_t*&) const;

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

		//span
		void span(sparse&, uint32_t, uint32_t, uint32_t, uint32_t) const;
		void span_data(sparse&, uint32_t, uint32_t, uint32_t, uint32_t) const;
		void span_pattern(sparse&, uint32_t, uint32_t, uint32_t, uint32_t) const;

	private:
		//data
		void cleanup(void);

		//print
		void print_dense(void) const;
		void print_sparse(void) const;

		//data
		bool m_own;
		double* m_data_ptr;
		const double* m_data_ref;

		uint32_t m_rows;
		uint32_t m_cols;
		int32_t* m_rows_map_ptr;
		int32_t* m_cols_map_ptr;
		const int32_t* m_rows_map_ref;
		const int32_t* m_cols_map_ref;
	};
}