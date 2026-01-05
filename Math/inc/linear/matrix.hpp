#pragma once

//std
#include <cstdint>
#include <initializer_list>

//math
#include "Math/Math/inc/linear/span.hpp"

//defines
#ifndef MATRIX_STATIC_SIZE
#define MATRIX_STATIC_SIZE 324
#endif

namespace math
{
	class mat3;
	class vector;
	enum class mode : uint32_t
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
		matrix(const double*, uint32_t, uint32_t);
		matrix(uint32_t, uint32_t, mode = mode::null);
		matrix(double*, uint32_t, uint32_t, mode = mode::null);
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

		double& operator[](uint32_t);
		double& operator()(uint32_t);
		double& operator()(uint32_t, uint32_t);

		const double& operator[](uint32_t) const;
		const double& operator()(uint32_t) const;
		const double& operator()(uint32_t, uint32_t) const;

		//data
		uint32_t rows(void) const;
		uint32_t cols(void) const;

		double* data(void);
		const double* data(void) const;

		//size
		matrix& resize(uint32_t, uint32_t);

		//bounds
		double min(bool = false, uint32_t* = nullptr) const;
		double max(bool = false, uint32_t* = nullptr) const;

		//util
		matrix& eye(void);
		matrix& zeros(void);
		matrix& swap_rows(uint32_t, uint32_t);
		matrix& swap_cols(uint32_t, uint32_t);
		matrix& randu(double = -1, double = +1);

		//info
		void print(const char* = "", double = 0) const;

		//linear
		double norm(void) const;
		double trace(void) const;
		double determinant(void) const;
		double bilinear(const double*) const;
		double bilinear(const vector&) const;
		double bilinear(const double*, const double*) const;
		double bilinear(const vector&, const vector&) const;

		matrix transpose(void) const;
		matrix inverse(bool* = nullptr) const;
		bool solve(matrix&, const matrix&) const;

		bool solve_decompose(uint32_t*);
		bool solve_substitute(const uint32_t*, matrix&);
		bool solve_substitute(const uint32_t*, double*, uint32_t);
		bool solve_substitute(const uint32_t*, const matrix&, matrix&);
		bool solve_substitute(const uint32_t*, const double*, double*, uint32_t);

		bool symmetric(double = 1e-5) const;

		mat3 span3(uint32_t, uint32_t) const;
		math::span span(uint32_t, uint32_t, uint32_t = 3, uint32_t = 3);

		//svd
		bool svd(matrix&, matrix&, vector&) const;

		//eigen
		bool eigen(vector&, vector&, matrix&) const;

		bool eigen_sym(vector&, matrix&) const;
		bool eigen_sym(vector&, matrix&, const matrix&) const;

		//stats
		double mean(void) const;
		double variance(void) const;

		//static
		static matrix eye(uint32_t, uint32_t);
		static matrix zeros(uint32_t, uint32_t);

		//friends
		friend matrix operator*(double, const matrix&);
		friend matrix operator*(const matrix&, const math::span&);
		friend matrix operator*(const math::span&, const matrix&);

	protected:
		//data
		bool m_own;
		uint32_t m_rows;
		uint32_t m_cols;
		double* m_data_ptr;
		const double* m_data_ref;
		double m_data_mem[MATRIX_STATIC_SIZE];

		//friends
		friend class span;
	};
}