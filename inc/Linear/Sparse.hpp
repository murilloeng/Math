#pragma once

//std
#include <cstdint>

namespace math
{
	class Vector;
	class Matrix;
}

namespace math
{
	class Sparse
	{
	public:
		//constructors
		Sparse(uint32_t, uint32_t);
		Sparse(double*, int32_t*, int32_t*, uint32_t, uint32_t);
		Sparse(double*, const int32_t*, const int32_t*, uint32_t, uint32_t);
		Sparse(const double*, const int32_t*, const int32_t*, uint32_t, uint32_t);

		//destructor
		~Sparse(void);

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

		//solve
		bool solve(Vector&, const Vector&) const;

		//print
		void print(const char* = "", bool = true) const;

		//operators
		Vector operator*(const Vector&) const;
		double& operator()(uint32_t, uint32_t);
		const double operator()(uint32_t, uint32_t) const;

		//convert
		Matrix convert(void) const;

		//Span
		void Span(Sparse&, uint32_t, uint32_t) const;
		void span_data(Sparse&, uint32_t, uint32_t) const;
		void span_pattern(Sparse&, uint32_t, uint32_t) const;

	private:
		//data
		void cleanup(void);

		//print
		void print_dense(void) const;
		void print_sparse(void) const;

		//Span
		void span_count(Sparse&, uint32_t, uint32_t, bool) const;
		void span_check(Sparse&, uint32_t, uint32_t, bool) const;

		//data
		bool m_own;
		double* m_data_ptr;
		const double* m_data_ref;

		mutable void* m_numeric;
		mutable void* m_symbolic;

		uint32_t m_rows;
		uint32_t m_cols;
		int32_t* m_rows_map_ptr;
		int32_t* m_cols_map_ptr;
		const int32_t* m_rows_map_ref;
		const int32_t* m_cols_map_ref;
	};
}