#pragma once

//std
#include <initializer_list>

//mat
#include "inc/linear/span.hpp"

//defines
#ifndef MATRIX_STATIC_SIZE
#define MATRIX_STATIC_SIZE 324
#endif

namespace math
{
	class mat3;
	class vector;
	enum class mode : unsigned
	{
		eye,
		null,
		zeros
	};
}

namespace math
{
	class matrix
	{
	public:
		//constructors
		matrix(void);
		matrix(const matrix&);
		matrix(const double*, unsigned, unsigned);
		matrix(unsigned, unsigned, mode = mode::null);
		matrix(double*, unsigned, unsigned, mode = mode::null);
		matrix(std::initializer_list<std::initializer_list<double>>, bool = false);

		//destructor
		virtual ~matrix(void);

		//serialization
		void load(const char*);
		void save(const char*) const;

		//operators
		matrix operator+(void) const;
		matrix operator-(void) const;
		matrix operator/(double) const;

		matrix operator+(const matrix&) const;
		matrix operator-(const matrix&) const;
		vector operator*(const vector&) const;
		matrix operator*(const matrix&) const;

		matrix& operator=(double);
		matrix& operator=(const double*);
		matrix& operator=(const matrix&);
		matrix& operator=(std::initializer_list<double>);
		matrix& operator=(std::initializer_list<std::initializer_list<double>>);

		matrix& operator+=(double);
		matrix& operator-=(double);
		matrix& operator*=(double);
		matrix& operator/=(double);

		matrix& operator+=(const double*);
		matrix& operator-=(const double*);
		matrix& operator+=(const matrix&);
		matrix& operator-=(const matrix&);

		double& operator[](unsigned);
		double& operator()(unsigned);
		double& operator()(unsigned, unsigned);

		const double& operator[](unsigned) const;
		const double& operator()(unsigned) const;
		const double& operator()(unsigned, unsigned) const;

		//data
		unsigned rows(void) const;
		unsigned cols(void) const;

		double* mem(void);
		const double* mem(void) const;

		//size
		matrix& resize(unsigned, unsigned);

		//bounds
		double min(bool = false, unsigned* = nullptr) const;
		double max(bool = false, unsigned* = nullptr) const;

		//util
		matrix& eye(void);
		matrix& zeros(void);
		matrix& swap_rows(unsigned, unsigned);
		matrix& swap_cols(unsigned, unsigned);
		matrix& randu(double = -1, double = +1);

		//info
		void print(const char* = "", double = 0) const;

		//linear
		double norm(void) const;
		double trace(void) const;
		double determinant(void) const;

		matrix inverse(void) const;
		matrix transpose(void) const;
		void solve(matrix&, const matrix&) const;

		bool symmetric(double = 1e-5) const;

		mat3 span3(unsigned, unsigned) const;
		math::span span(unsigned, unsigned, unsigned = 3, unsigned = 3);

		//eigen
		bool eigen_sym(vector&, matrix&) const;
		bool eigen_sym(vector&, matrix&, const matrix&) const;

		//stats
		double mean(void) const;
		double variance(void) const;

		//static
		static matrix eye(unsigned, unsigned);
		static matrix zeros(unsigned, unsigned);

		//friends
		friend matrix operator*(double, const matrix&);
		friend matrix operator*(const matrix&, const math::span&);
		friend matrix operator*(const math::span&, const matrix&);

	protected:
		//data
		bool m_own;
		double* m_ptr;
		unsigned m_rows;
		unsigned m_cols;
		const double* m_ref;
		double m_mem[MATRIX_STATIC_SIZE];

		//friends
		friend class span;
	};
}