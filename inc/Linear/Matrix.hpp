#pragma once

//std
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <concepts>

//Math
#include "Math/inc/Linear/Modes.hpp"

extern "C"
{
	void dgetrf_(const uint32_t*, const uint32_t*, double*, const uint32_t*, uint32_t*, int32_t*);
	void dgetri_(const uint32_t*, double*, const uint32_t*, const uint32_t*, double*, const int32_t*, int32_t*);
}

namespace math
{
	template<uint32_t r, uint32_t c, uint32_t m, uint32_t n> requires (m <= r && n <= c) class SpanFixed;
}

namespace math
{
	template<uint32_t r, uint32_t c> class Matrix
	{
	public:
		//constructor
		Matrix(Mode = Mode::none);
		explicit Matrix(const double*);
		template<class... T> Matrix(T...) requires (sizeof...(T) == r * c && (std::convertible_to<T, double> && ...));
		template<uint32_t... k> Matrix(const double (&... rows)[k]) requires (sizeof...(k) == r && ((k == c) && ...));

		//destructor
		~Matrix(void) = default;

		//write
		void write(double*) const;

		//data
		double* data(void);
		const double* data(void) const;

		//operators
		Matrix operator+(void) const;
		Matrix operator-(void) const;

		Matrix& operator+=(double);
		Matrix& operator-=(double);
		Matrix& operator*=(double);
		Matrix& operator/=(double);
		Matrix& operator+=(const Matrix&);
		Matrix& operator-=(const Matrix&);

		Matrix operator/(double) const;
		Matrix operator+(const Matrix&) const;
		Matrix operator-(const Matrix&) const;

		double& operator[](uint32_t);
		double& operator()(uint32_t, uint32_t);
		const double& operator[](uint32_t) const;
		const double& operator()(uint32_t, uint32_t) const;

		template<uint32_t d> Matrix<r, d> operator*(const Matrix<c, d>&) const;

		//linear
		Matrix& eye(void);
		Matrix& zeros(void);
		Matrix& randu(double = -1, double = +1);
		Matrix& normalize(void);

		Matrix unit(void) const;
		Matrix<c, r> transpose(void) const;
		Matrix inverse(bool* = nullptr) const requires (r == c);

		double norm(void) const;
		double inner(void) const requires (c == 1);
		double inner(const Matrix&) const requires (c == 1);
		double bilinear(const Matrix<c, 1>&) const requires (r == c);
		double bilinear(const Matrix<r, 1>&, const Matrix<c, 1>&) const;

		Matrix<r, r> spin(void) const requires (r == 3 && c == 1);
		double cross(const Matrix&) const requires (r == 2 && c == 1);
		Matrix cross(const Matrix&) const requires (r == 3 && c == 1);
		void triad(Matrix&, Matrix&) const requires (r == 3 && c == 1);

		Matrix<r, r> outer(void) const requires (c == 1);
		template<uint32_t a> Matrix<r, a> outer(const Matrix<a, 1>&) const requires (c == 1);

		//span
		template<uint32_t m, uint32_t n> SpanFixed<r, c, m, n> span(uint32_t, uint32_t) requires (m <= r && n <= c);

		//friends
		template<uint32_t a, uint32_t b> friend Matrix<a, b> operator*(double, const Matrix<a, b>&);

	protected:
		//data
		double m_data[r * c];

		//friends
		template<uint32_t, uint32_t> friend class Matrix;
	};

	//alias
	using Vec1 = Matrix<1, 1>;
	using Vec2 = Matrix<2, 1>;
	using Vec3 = Matrix<3, 1>;
	using Vec4 = Matrix<4, 1>;
	using Vec6 = Matrix<6, 1>;
	using Mat1 = Matrix<1, 1>;
	using Mat2 = Matrix<2, 2>;
	using Mat3 = Matrix<3, 3>;
	using Mat4 = Matrix<4, 4>;
	using Mat6 = Matrix<6, 6>;
	template<uint32_t n> using Vector = Matrix<n, 1>;

	//constructors
	template<uint32_t r, uint32_t c> Matrix<r, c>::Matrix(Mode mode)
	{
		if(mode == Mode::eye) eye();
		if(mode == Mode::zeros) zeros();
	}
	template<uint32_t r, uint32_t c> Matrix<r, c>::Matrix(const double* data)
	{
		memcpy(m_data, data, r * c * sizeof(double));
	}
	template<uint32_t r, uint32_t c> template<class... T> Matrix<r, c>::Matrix(T... data) requires (sizeof...(T) == r * c && (std::convertible_to<T, double> && ...))
	{
		const double list[] = {double(data)...};
		for(uint32_t i = 0; i < r; i++)
		{
			for(uint32_t j = 0; j < c; j++)
			{
				m_data[i + r * j] = list[c * i + j];
			}
		}
	}
	template<uint32_t r, uint32_t c> template<uint32_t... k> Matrix<r, c>::Matrix(const double (&... rows)[k]) requires (sizeof...(k) == r && ((k == c) && ...))
	{
		const double* list[] = {rows...};
		for(uint32_t i = 0; i < r; i++)
		{
			for(uint32_t j = 0; j < c; j++)
			{
				m_data[i + r * j] = list[i][j];
			}
		}
	}

	//write
	template<uint32_t r, uint32_t c> void Matrix<r, c>::write(double* data) const
	{
		memcpy(data, m_data, r * c * sizeof(double));
	}

	//data
	template<uint32_t r, uint32_t c> double* Matrix<r, c>::data(void)
	{
		return m_data;
	}
	template<uint32_t r, uint32_t c> const double* Matrix<r, c>::data(void) const
	{
		return m_data;
	}

	//operators
	template<uint32_t r, uint32_t c> Matrix<r, c> Matrix<r, c>::operator+(void) const
	{
		return *this;
	}
	template<uint32_t r, uint32_t c> Matrix<r, c> Matrix<r, c>::operator-(void) const
	{
		Matrix R;
		for(uint32_t i = 0; i < r * c; i++) R.m_data[i] = -m_data[i];
		return R;
	}

	template<uint32_t r, uint32_t c> Matrix<r, c>& Matrix<r, c>::operator+=(double s)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data[i] += s;
		return *this;
	}
	template<uint32_t r, uint32_t c> Matrix<r, c>& Matrix<r, c>::operator-=(double s)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data[i] -= s;
		return *this;
	}
	template<uint32_t r, uint32_t c> Matrix<r, c>& Matrix<r, c>::operator*=(double s)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data[i] *= s;
		return *this;
	}
	template<uint32_t r, uint32_t c> Matrix<r, c>& Matrix<r, c>::operator/=(double s)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data[i] /= s;
		return *this;
	}
	template<uint32_t r, uint32_t c> Matrix<r, c>& Matrix<r, c>::operator+=(const Matrix& M)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data[i] += M.m_data[i];
		return *this;
	}
	template<uint32_t r, uint32_t c> Matrix<r, c>& Matrix<r, c>::operator-=(const Matrix& M)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data[i] -= M.m_data[i];
		return *this;
	}

	template<uint32_t r, uint32_t c>  Matrix<r, c> Matrix<r, c>::operator/(double s) const
	{
		Matrix R;
		for(uint32_t i = 0; i < r * c; i++) R.m_data[i] = m_data[i] / s;
		return R;
	}
	template<uint32_t r, uint32_t c>  Matrix<r, c> Matrix<r, c>::operator+(const Matrix& M) const
	{
		Matrix R;
		for(uint32_t i = 0; i < r * c; i++) R.m_data[i] = m_data[i] + M.m_data[i];
		return R;
	}
	template<uint32_t r, uint32_t c> Matrix<r, c> Matrix<r, c>::operator-(const Matrix& M) const
	{
		Matrix R;
		for(uint32_t i = 0; i < r * c; i++) R.m_data[i] = m_data[i] - M.m_data[i];
		return R;
	}

	template<uint32_t r, uint32_t c> double& Matrix<r, c>::operator[](uint32_t index)
	{
		return m_data[index];
	}
	template<uint32_t r, uint32_t c> double& Matrix<r, c>::operator()(uint32_t index_1, uint32_t index_2)
	{
		return m_data[index_1 + r * index_2];
	}
	template<uint32_t r, uint32_t c> const double& Matrix<r, c>::operator[](uint32_t index) const
	{
		return m_data[index];
	}
	template<uint32_t r, uint32_t c> const double& Matrix<r, c>::operator()(uint32_t index_1, uint32_t index_2) const
	{
		return m_data[index_1 + r * index_2];
	}

	template<uint32_t r, uint32_t c> template<uint32_t d> Matrix<r, d> Matrix<r, c>::operator*(const Matrix<c, d>& M) const
	{
		double v[r];
		Matrix<r, d> R;
		for(uint32_t j = 0; j < d; j++)
		{
			for(uint32_t i = 0; i < r; i++) v[i] = 0;
			for(uint32_t k = 0; k < c; k++)
			{
				const double s = M.m_data[k + c * j];
				for(uint32_t i = 0; i < r; i++) v[i] += m_data[i + r * k] * s;
			}
			for(uint32_t i = 0; i < r; i++) R.m_data[i + r * j] = v[i];
		}
		return R;
	}

	//linear
	template<uint32_t r, uint32_t c> Matrix<r, c>& Matrix<r, c>::eye(void)
	{
		for(uint32_t i = 0; i < r; i++)
		{
			for(uint32_t j = 0; j < c; j++)
			{
				m_data[i + r * j] = i == j;
			}
		}
		return *this;
	}
	template<uint32_t r, uint32_t c> Matrix<r, c>& Matrix<r, c>::zeros(void)
	{
		memset(m_data, 0, r * c * sizeof(double));
		return *this;
	}
	template<uint32_t r, uint32_t c> Matrix<r, c>& Matrix<r, c>::randu(double a, double b)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data[i] = a + (b - a) * rand() / RAND_MAX;
		return *this;
	}
	template<uint32_t r, uint32_t c> Matrix<r, c>& Matrix<r, c>::normalize(void)
	{
		return *this /= norm();
	}

	template<uint32_t r, uint32_t c> Matrix<r, c> Matrix<r, c>::unit(void) const
	{
		return *this / norm();
	}
	template<uint32_t r, uint32_t c> Matrix<c, r> Matrix<r, c>::transpose(void) const
	{
		Matrix<c, r> R;
		for(uint32_t i = 0; i < c; i++)
		{
			for(uint32_t j = 0; j < r; j++)
			{
				R.m_data[i + c * j] = m_data[j + r * i];
			}
		}
		return R;
	}
	template<uint32_t r, uint32_t c> Matrix<r, c> Matrix<r, c>::inverse(bool* test) const requires (r == c)
	{
		//query
		double query;
		const uint32_t n = r;
		int32_t info, lwork = -1;
		Matrix<r, c> A(*this);
		uint32_t* ipiv = new uint32_t[n];
		dgetrf_(&n, &n, A.m_data, &n, ipiv, &info);
		dgetri_(&n, A.m_data, &n, ipiv, &query, &lwork, &info);
		//inverse
		lwork = int32_t(query);
		double* work = new double[lwork];
		dgetri_(&n, A.m_data, &n, ipiv, work, &lwork, &info);
		//delete
		delete[] ipiv;
		delete[] work;
		if(test) *test = info == 0;
		//return
		return A;
	}

	template<uint32_t r, uint32_t c> double Matrix<r, c>::norm(void) const
	{
		double s = 0;
		for(uint32_t i = 0; i < r * c; i++)
		{
			s += m_data[i] * m_data[i];
		}
		return sqrt(s);
	}
	template<uint32_t r, uint32_t c> double Matrix<r, c>::inner(void) const requires (c == 1)
	{
		return inner(*this);
	}
	template<uint32_t r, uint32_t c> double Matrix<r, c>::inner(const Matrix& v) const requires (c == 1)
	{
		double s = 0;
		for(uint32_t i = 0; i < r; i++)
		{
			s += m_data[i] * v.m_data[i];
		}
		return s;
	}
	template<uint32_t r, uint32_t c> double Matrix<r, c>::bilinear(const Matrix<c, 1>& v) const requires (r == c)
	{
		return bilinear(v, v);
	}
	template<uint32_t r, uint32_t c> double Matrix<r, c>::bilinear(const Matrix<r, 1>& v1, const Matrix<c, 1>& v2) const
	{
		double s = 0;
		for(uint32_t i = 0; i < r; i++)
		{
			for(uint32_t j = 0; j < c; j++)
			{
				s += v1.m_data[i] * m_data[i + r * j] * v2.m_data[j];
			}
		}
		return s;
	}

	template<uint32_t r, uint32_t c> Matrix<r, r> Matrix<r, c>::spin(void) const requires (r == 3 && c == 1)
	{
		Matrix<r, r> R;
		R.m_data[7] = -m_data[0];
		R.m_data[2] = -m_data[1];
		R.m_data[3] = -m_data[2];
		R.m_data[5] = +m_data[0];
		R.m_data[6] = +m_data[1];
		R.m_data[1] = +m_data[2];
		R.m_data[0] = R.m_data[4] = R.m_data[8] = 0;
		return R;
	}
	template<uint32_t r, uint32_t c> double Matrix<r, c>::cross(const Matrix& v) const requires (r == 2 && c == 1)
	{
		return m_data[0] * v.m_data[1] - m_data[1] * v.m_data[0];
	}
	template<uint32_t r, uint32_t c> Matrix<r, c> Matrix<r, c>::cross(const Matrix& v) const requires (r == 3 && c == 1)
	{
		Matrix<r, c> w;
		w.m_data[0] = m_data[1] * v.m_data[2] - m_data[2] * v.m_data[1];
		w.m_data[1] = m_data[2] * v.m_data[0] - m_data[0] * v.m_data[2];
		w.m_data[2] = m_data[0] * v.m_data[1] - m_data[1] * v.m_data[0];
		return w;
	}
	template<uint32_t r, uint32_t c> void Matrix<r, c>::triad(Matrix& v2, Matrix& v3) const requires (r == 3 && c == 1)
	{
		//index
		uint32_t index = 0;
		const Vec3 v1 = unit();
		if(fabs(v1.m_data[1]) < fabs(v1.m_data[index])) index = 1;
		if(fabs(v1.m_data[2]) < fabs(v1.m_data[index])) index = 2;
		//normal
		v2.m_data[index] = 0;
		v2.m_data[(index + 1) % 3] = +v1.m_data[(index + 2) % 3];
		v2.m_data[(index + 2) % 3] = -v1.m_data[(index + 1) % 3];
		//triad
		v2.normalize();
		v3 = v1.cross(v2);
	}

	template<uint32_t r, uint32_t c> Matrix<r, r> Matrix<r, c>::outer(void) const requires (c == 1)
	{
		return outer(*this);
	}
	template<uint32_t r, uint32_t c> template<uint32_t a> Matrix<r, a> Matrix<r, c>::outer(const Matrix<a, 1>& v) const requires (c == 1)
	{
		Matrix<r, a> R;
		for(uint32_t i = 0; i < r; i++)
		{
			for(uint32_t j = 0; j < a; j++)
			{
				R.m_data[i + r * j] = m_data[i] * v.m_data[j];
			}
		}
		return R;
	}

	//span
	template<uint32_t r, uint32_t c> template<uint32_t m, uint32_t n> SpanFixed<r, c, m, n> Matrix<r, c>::span(uint32_t row, uint32_t col) requires (m <= r && n <= c)
	{
		return SpanFixed<r, c, m, n>(*this, row, col);
	}

	//friends
	template<uint32_t r, uint32_t c> Matrix<r, c> operator*(double s, const Matrix<r, c>& M)
	{
		Matrix<r, c> R;
		for(uint32_t i = 0; i < r * c; i++) R.m_data[i] = s * M.m_data[i];
		return R;
	}
}

//Math
#include "Math/inc/Linear/SpanFixed.hpp"