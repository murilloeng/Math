#pragma once

//std
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <initializer_list>

//Math
#include "Math/inc/Linear/Modes.hpp"
#include "Math/inc/Linear/Matrix.hpp"

//defines
#ifndef MATRIX_STATIC_SIZE
#define MATRIX_STATIC_SIZE 324
#endif

namespace math
{
	class SpanN;
	class VectorN;
}

namespace math
{
	class MatrixN
	{
	public:
		//constructors
		MatrixN(const MatrixN&);
		MatrixN(uint32_t, uint32_t);
		MatrixN(uint32_t, uint32_t, Mode);
		MatrixN(double*, uint32_t, uint32_t);
		MatrixN(const double*, uint32_t, uint32_t);
		MatrixN(std::initializer_list<std::initializer_list<double>>);

		//destructor
		virtual ~MatrixN(void);

		//operators
		MatrixN operator+(void) const;
		MatrixN operator-(void) const;
		MatrixN operator/(double) const;

		MatrixN operator+(const MatrixN&) const;
		MatrixN operator-(const MatrixN&) const;
		VectorN operator*(const VectorN&) const;
		MatrixN operator*(const MatrixN&) const;

		MatrixN& operator=(double);
		MatrixN& operator=(const double*);
		MatrixN& operator=(const MatrixN&);
		MatrixN& operator=(std::initializer_list<double>);
		MatrixN& operator=(std::initializer_list<std::initializer_list<double>>);
		template<uint32_t r, uint32_t c> MatrixN& operator=(const Matrix<r, c>&);

		MatrixN& operator+=(double);
		MatrixN& operator-=(double);
		MatrixN& operator*=(double);
		MatrixN& operator/=(double);

		MatrixN& operator+=(const MatrixN&);
		MatrixN& operator-=(const MatrixN&);
		MatrixN& operator+=(const double*);
		MatrixN& operator-=(const double*);

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

		//bounds
		double min(bool = false, uint32_t* = nullptr) const;
		double max(bool = false, uint32_t* = nullptr) const;

		//utility
		MatrixN& eye(void);
		MatrixN& zeros(void);
		MatrixN& randu(double = -1, double = +1);

		//print
		void print(const char* = "", double = 0) const;

		//linear
		double norm(void) const;
		double trace(void) const;
		double determinant(void) const;
		double bilinear(const double*) const;
		double bilinear(const VectorN&) const;
		double bilinear(const double*, const double*) const;
		double bilinear(const VectorN&, const VectorN&) const;

		MatrixN transpose(void) const;
		MatrixN inverse(bool* = nullptr) const;
		bool solve(MatrixN&, const MatrixN&) const;
		bool solve(VectorN&, const VectorN&) const;
		bool solve(double*, const double*, uint32_t = 1) const;

		bool solve_decompose(uint32_t*);
		bool solve_substitute(const uint32_t*, MatrixN&);
		bool solve_substitute(const uint32_t*, double*, uint32_t);
		bool solve_substitute(const uint32_t*, const MatrixN&, MatrixN&);
		bool solve_substitute(const uint32_t*, const double*, double*, uint32_t);

		SpanN span(uint32_t, uint32_t, uint32_t = 3, uint32_t = 3);

		//stats
		double mean(void) const;
		double variance(void) const;

		//static
		static MatrixN eye(uint32_t, uint32_t);
		static MatrixN zeros(uint32_t, uint32_t);

		//friends
		friend MatrixN operator*(double, const MatrixN&);
		friend MatrixN operator*(const MatrixN&, const SpanN&);
		friend MatrixN operator*(const SpanN&, const MatrixN&);

	protected:
		//data
		bool m_own;
		uint32_t m_rows;
		uint32_t m_cols;
		double* m_data_ptr;
		const double* m_data_ref;
		double m_data_mem[MATRIX_STATIC_SIZE];

		//friends
		friend class math::SpanN;
		friend class math::VectorN;
	};
}