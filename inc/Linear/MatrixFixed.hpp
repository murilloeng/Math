#pragma once

//std
#include <cstdint>
#include <cstring>
#include <concepts>

namespace math
{
	template<uint32_t r, uint32_t c> class MatrixFixed
	{
	public:
		//constructor
		MatrixFixed(void) = default;
		explicit MatrixFixed(const double*);
		template<class... T> MatrixFixed(T...) requires (sizeof...(T) == r * c && (std::convertible_to<T, double> && ...));
		template<uint32_t... k> MatrixFixed(const double (&... rows)[k]) requires (sizeof...(k) == r && ((k == c) && ...));

		//destructor
		~MatrixFixed(void) = default;

		//data
		double* data(void);
		const double* data(void) const;

		//operators
		MatrixFixed operator+(void) const;
		MatrixFixed operator-(void) const;

		MatrixFixed& operator+=(double);
		MatrixFixed& operator-=(double);
		MatrixFixed& operator*=(double);
		MatrixFixed& operator/=(double);
		MatrixFixed& operator+=(const MatrixFixed&);
		MatrixFixed& operator-=(const MatrixFixed&);

		MatrixFixed operator/(double) const;
		MatrixFixed operator+(const MatrixFixed&) const;
		MatrixFixed operator-(const MatrixFixed&) const;

		double& operator[](uint32_t);
		double& operator()(uint32_t, uint32_t);
		const double& operator[](uint32_t) const;
		const double& operator()(uint32_t, uint32_t) const;

		template<uint32_t d> MatrixFixed<r, d> operator*(const MatrixFixed<c, d>&) const;

		//linear
		MatrixFixed& eye(void);
		MatrixFixed& zeros(void);
		MatrixFixed<c, r> transpose(void) const;
		double bilinear(const MatrixFixed<c, 1>&) const requires (r == c);
		double bilinear(const MatrixFixed<r, 1>&, const MatrixFixed<c, 1>&) const;

		//friends
		template<uint32_t a, uint32_t b> friend MatrixFixed<a, b> operator*(double, const MatrixFixed<a, b>&);

	protected:
		//data
		double m_data[r * c] = {};

		//friends
		template<uint32_t, uint32_t> friend class MatrixFixed;
	};

	//constructors
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>::MatrixFixed(const double* data)
	{
		memcpy(m_data, data, r * c * sizeof(double));
	}
	template<uint32_t r, uint32_t c> template<class... T> MatrixFixed<r, c>::MatrixFixed(T... data) requires (sizeof...(T) == r * c && (std::convertible_to<T, double> && ...))
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
	template<uint32_t r, uint32_t c> template<uint32_t... k> MatrixFixed<r, c>::MatrixFixed(const double (&... rows)[k]) requires (sizeof...(k) == r && ((k == c) && ...))
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

	//data
	template<uint32_t r, uint32_t c>  double* MatrixFixed<r, c>::data(void)
	{
		return m_data;
	}
	template<uint32_t r, uint32_t c> const double* MatrixFixed<r, c>::data(void) const
	{
		return m_data;
	}

	//operators
	template<uint32_t r, uint32_t c> MatrixFixed<r, c> MatrixFixed<r, c>::operator+(void) const
	{
		return *this;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c> MatrixFixed<r, c>::operator-(void) const
	{
		MatrixFixed R;
		for(uint32_t i = 0; i < r * c; i++) R.m_data[i] = -m_data[i];
		return R;
	}

	template<uint32_t r, uint32_t c> MatrixFixed<r, c>& MatrixFixed<r, c>::operator+=(double s)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data[i] += s;
		return *this;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>& MatrixFixed<r, c>::operator-=(double s)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data[i] -= s;
		return *this;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>& MatrixFixed<r, c>::operator*=(double s)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data[i] *= s;
		return *this;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>& MatrixFixed<r, c>::operator/=(double s)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data[i] /= s;
		return *this;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>& MatrixFixed<r, c>::operator+=(const MatrixFixed& M)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data[i] += M.m_data[i];
		return *this;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>& MatrixFixed<r, c>::operator-=(const MatrixFixed& M)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data[i] -= M.m_data[i];
		return *this;
	}

	template<uint32_t r, uint32_t c>  MatrixFixed<r, c> MatrixFixed<r, c>::operator/(double s) const
	{
		MatrixFixed R;
		for(uint32_t i = 0; i < r * c; i++) R.m_data[i] = m_data[i] / s;
		return R;
	}
	template<uint32_t r, uint32_t c>  MatrixFixed<r, c> MatrixFixed<r, c>::operator+(const MatrixFixed& M) const
	{
		MatrixFixed R;
		for(uint32_t i = 0; i < r * c; i++) R.m_data[i] = m_data[i] + M.m_data[i];
		return R;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c> MatrixFixed<r, c>::operator-(const MatrixFixed& M) const
	{
		MatrixFixed R;
		for(uint32_t i = 0; i < r * c; i++) R.m_data[i] = m_data[i] - M.m_data[i];
		return R;
	}

	template<uint32_t r, uint32_t c> double& MatrixFixed<r, c>::operator[](uint32_t index)
	{
		return m_data[index];
	}
	template<uint32_t r, uint32_t c> double& MatrixFixed<r, c>::operator()(uint32_t index_1, uint32_t index_2)
	{
		return m_data[index_1 + r * index_2];
	}
	template<uint32_t r, uint32_t c> const double& MatrixFixed<r, c>::operator[](uint32_t index) const
	{
		return m_data[index];
	}
	template<uint32_t r, uint32_t c> const double& MatrixFixed<r, c>::operator()(uint32_t index_1, uint32_t index_2) const
	{
		return m_data[index_1 + r * index_2];
	}

	template<uint32_t r, uint32_t c> template<uint32_t d> MatrixFixed<r, d> MatrixFixed<r, c>::operator*(const MatrixFixed<c, d>& M) const
	{
		double v[r];
		MatrixFixed<r, d> R;
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
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>& MatrixFixed<r, c>::eye(void)
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
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>& MatrixFixed<r, c>::zeros(void)
	{
		memset(m_data, 0, r * c * sizeof(double));
		return *this;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<c, r> MatrixFixed<r, c>::transpose(void) const
	{
		MatrixFixed<c, r> R;
		for(uint32_t i = 0; i < c; i++)
		{
			for(uint32_t j = 0; j < r; j++)
			{
				R.m_data[i + c * j] = m_data[j + r * i];
			}
		}
		return R;
	}
	template<uint32_t r, uint32_t c> double MatrixFixed<r, c>::bilinear(const MatrixFixed<c, 1>& v) const requires (r == c)
	{
		return bilinear(v, v);
	}
	template<uint32_t r, uint32_t c> double MatrixFixed<r, c>::bilinear(const MatrixFixed<r, 1>& v1, const MatrixFixed<c, 1>& v2) const
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

	//friends
	template<uint32_t r, uint32_t c> MatrixFixed<r, c> operator*(double s, const MatrixFixed<r, c>& M)
	{
		MatrixFixed<r, c> R;
		for(uint32_t i = 0; i < r * c; i++) R.m_data[i] = s * M.m_data[i];
		return R;
	}
}