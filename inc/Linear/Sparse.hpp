#pragma once

//std
#include <cstdint>

namespace math
{
	class VectorX;
	class MatrixX;
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
		void zeros(void);
		double norm(void) const;
		double trace(void) const;
		double bilinear(const double*) const;
		double bilinear(const VectorX&) const;
		double bilinear(const double*, const double*) const;
		double bilinear(const VectorX&, const VectorX&) const;

		//solve
		bool solve(VectorX&, const VectorX&) const;
		bool solve(MatrixX&, const MatrixX&) const;
		bool solve(double*, const double*, uint32_t) const;

		void release(void) const;
		bool decompose(void) const;
		bool substitute(VectorX&, const VectorX&) const;
		bool substitute(MatrixX&, const MatrixX&) const;
		bool substitute(double*, const double*, uint32_t) const;

		//operators
		VectorX operator*(const VectorX&) const;
		double& operator()(uint32_t, uint32_t);
		const double& operator()(uint32_t, uint32_t) const;

		//convert
		MatrixX convert(void) const;

		//product
		void product(double*, const double*) const;

		//span
		void span(Sparse&, uint32_t, uint32_t) const;
		void span_data(Sparse&, uint32_t, uint32_t) const;
		void span_pattern(Sparse&, uint32_t, uint32_t) const;

		//print
		void print(const char* = "", bool = true, double = 0) const;

	private:
		//data
		void cleanup(void);

		//print
		void print_dense(double) const;
		void print_sparse(double) const;

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