#pragma once

//std
#include <cstdint>
#include <initializer_list>

//Math
#include "Math/inc/Linear/Modes.hpp"

//defines
#ifndef MATRIX_STATIC_SIZE
#define MATRIX_STATIC_SIZE 324
#endif

namespace math
{
	class SpanX;
	class VectorX;
}

namespace math
{
	class MatrixX
	{
	public:
		//constructors
		MatrixX(const MatrixX&);
		MatrixX(uint32_t, uint32_t);
		MatrixX(uint32_t, uint32_t, Mode);
		MatrixX(double*, uint32_t, uint32_t);
		MatrixX(const double*, uint32_t, uint32_t);
		MatrixX(std::initializer_list<std::initializer_list<double>>);

		//destructor
		virtual ~MatrixX(void);

		//operators
		MatrixX operator+(void) const;
		MatrixX operator-(void) const;
		MatrixX operator/(double) const;

		MatrixX operator+(const MatrixX&) const;
		MatrixX operator-(const MatrixX&) const;
		VectorX operator*(const VectorX&) const;
		MatrixX operator*(const MatrixX&) const;

		MatrixX& operator=(double);
		MatrixX& operator=(const double*);
		MatrixX& operator=(const MatrixX&);
		MatrixX& operator=(std::initializer_list<double>);
		MatrixX& operator=(std::initializer_list<std::initializer_list<double>>);

		MatrixX& operator+=(double);
		MatrixX& operator-=(double);
		MatrixX& operator*=(double);
		MatrixX& operator/=(double);

		MatrixX& operator+=(const MatrixX&);
		MatrixX& operator-=(const MatrixX&);
		MatrixX& operator+=(const double*);
		MatrixX& operator-=(const double*);

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
		MatrixX& eye(void);
		MatrixX& zeros(void);
		MatrixX& randu(double = -1, double = +1);

		//print
		void print(const char* = "", double = 0) const;

		//linear
		double norm(void) const;
		double trace(void) const;
		double determinant(void) const;
		double bilinear(const double*) const;
		double bilinear(const VectorX&) const;
		double bilinear(const double*, const double*) const;
		double bilinear(const VectorX&, const VectorX&) const;

		MatrixX transpose(void) const;
		MatrixX inverse(bool* = nullptr) const;
		bool solve(MatrixX&, const MatrixX&) const;
		bool solve(VectorX&, const VectorX&) const;
		bool solve(double*, const double*, uint32_t = 1) const;

		bool solve_decompose(uint32_t*);
		bool solve_substitute(const uint32_t*, MatrixX&);
		bool solve_substitute(const uint32_t*, double*, uint32_t);
		bool solve_substitute(const uint32_t*, const MatrixX&, MatrixX&);
		bool solve_substitute(const uint32_t*, const double*, double*, uint32_t);

		SpanX span(uint32_t, uint32_t, uint32_t = 3, uint32_t = 3);

		//stats
		double mean(void) const;
		double variance(void) const;

		//static
		static MatrixX eye(uint32_t, uint32_t);
		static MatrixX zeros(uint32_t, uint32_t);

		//friends
		friend MatrixX operator*(double, const MatrixX&);
		friend MatrixX operator*(const MatrixX&, const SpanX&);
		friend MatrixX operator*(const SpanX&, const MatrixX&);

	protected:
		//data
		bool m_own;
		uint32_t m_rows;
		uint32_t m_cols;
		double* m_data_ptr;
		const double* m_data_ref;
		double m_data_mem[MATRIX_STATIC_SIZE];

		//friends
		friend class math::SpanX;
		friend class math::VectorX;
	};
}