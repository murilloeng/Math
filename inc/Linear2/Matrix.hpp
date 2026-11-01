#pragma once

//std
#include <cstdint>
#include <initializer_list>

//Math
#include "Math/inc/Linear/Span.hpp"

namespace math
{
	namespace linear
	{
		class Vector;
	}
}

namespace math
{
	namespace linear
	{
		class Matrix
		{
		public:
			//constructors
			Matrix(const Matrix&);
			Matrix(uint32_t, uint32_t);
			Matrix(double*, uint32_t, uint32_t);
			Matrix(const double*, uint32_t, uint32_t);
			Matrix(std::initializer_list<std::initializer_list<double>>);
	
			//destructor
			virtual ~Matrix(void);
	
			//operators
			Matrix operator+(void) const;
			Matrix operator-(void) const;
			Matrix operator/(double) const;
	
			Matrix operator+(const Matrix&) const;
			Matrix operator-(const Matrix&) const;
			Vector operator*(const Vector&) const;
			Matrix operator*(const Matrix&) const;
	
			Matrix& operator=(double);
			Matrix& operator=(const double*);
			Matrix& operator=(const Matrix&);
			Matrix& operator=(std::initializer_list<double>);
			Matrix& operator=(std::initializer_list<std::initializer_list<double>>);
	
			Matrix& operator+=(double);
			Matrix& operator-=(double);
			Matrix& operator*=(double);
			Matrix& operator/=(double);
	
			Matrix& operator+=(const double*);
			Matrix& operator-=(const double*);
			Matrix& operator+=(const Matrix&);
			Matrix& operator-=(const Matrix&);
	
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
	
			//util
			Matrix& eye(void);
			Matrix& zeros(void);
			Matrix& randu(double = -1, double = +1);
	
			//info
			void print(const char* = "", double = 0) const;
	
			//linear
			double norm(void) const;
			double trace(void) const;
			double determinant(void) const;
			double bilinear(const double*) const;
			double bilinear(const Vector&) const;
			double bilinear(const double*, const double*) const;
			double bilinear(const Vector&, const Vector&) const;
	
			Matrix transpose(void) const;
			Matrix inverse(bool* = nullptr) const;
			bool solve(Matrix&, const Matrix&) const;
			bool solve(double*, const double*, uint32_t = 1) const;
	
			bool solve_decompose(uint32_t*);
			bool solve_substitute(const uint32_t*, Matrix&);
			bool solve_substitute(const uint32_t*, double*, uint32_t);
			bool solve_substitute(const uint32_t*, const Matrix&, Matrix&);
			bool solve_substitute(const uint32_t*, const double*, double*, uint32_t);
	
			math::Span span(uint32_t, uint32_t, uint32_t = 3, uint32_t = 3);
	
			//stats
			double mean(void) const;
			double variance(void) const;
	
			//static
			static Matrix eye(uint32_t, uint32_t);
			static Matrix zeros(uint32_t, uint32_t);
	
			//friends
			friend Matrix operator*(double, const Matrix&);
			friend Matrix operator*(const Matrix&, const math::Span&);
			friend Matrix operator*(const math::Span&, const Matrix&);
	
		protected:
			//data
			bool m_own;
			uint32_t m_rows;
			uint32_t m_cols;
			double* m_data_ptr;
			const double* m_data_ref;
	
			//friends
			friend class Span;
		};
	}
}