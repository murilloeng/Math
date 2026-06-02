//std
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Linear/Vector.hpp"

namespace math
{
	//constructors
	Vector::Vector(void) : Matrix(0, 1)
	{
		return;
	}
	Vector::Vector(uint32_t rows, mode init) : Matrix(rows, 1, init)
	{
		return;
	}
	Vector::Vector(const Matrix& m) : Matrix(m.rows() * m.cols(), 1)
	{
		memcpy(m_data_ptr, m.data(), m_rows * sizeof(double));
	}
	Vector::Vector(const double* ptr, uint32_t rows) : Matrix(ptr, rows, 1)
	{
		return;
	}
	Vector::Vector(std::initializer_list<double> list) : Matrix(list.size(), 1)
	{
		memcpy(m_data_ptr, std::data(list), list.size() * sizeof(double));
	}
	Vector::Vector(double* ptr, uint32_t rows, mode init) : Matrix(ptr, rows, 1, init)
	{
		return;
	}

	//destructor
	Vector::~Vector(void)
	{
		return;
	}

	//size
	Vector& Vector::resize(uint32_t rows)
	{
		Matrix::resize(rows, 1);
		return *this;
	}
	Vector& Vector::resize(uint32_t rows, uint32_t cols)
	{
		//check
		if(cols != 1)
		{
			throw std::runtime_error("Vector resize called with number of columns not equal to 1!");
		}
		//resize
		Matrix::resize(rows, 1);
		return *this;
	}

	//linear
	Vector& Vector::normalize(void)
	{
		return *this /= norm();
	}
	Vector Vector::unit(void) const
	{
		return *this / norm();
	}
	Matrix Vector::outer(void) const
	{
		return outer(*this);
	}
	Matrix Vector::outer(const Vector& v) const
	{
		Matrix m(m_rows, v.m_rows);
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < v.m_rows; j++)
			{
				m(i, j) = m_data_ref[i] * v.m_data_ref[j];
			}
		}
		return m;
	}
	double Vector::inner(const Vector& v) const
	{
		double s = 0;
		for(uint32_t i = 0; i < m_rows; i++)
		{
			s += m_data_ref[i] * v.m_data_ref[i];
		}
		return s;
	}
	double Vector::inner(const double* v) const
	{
		double s = 0;
		for(uint32_t i = 0; i < m_rows; i++)
		{
			s += m_data_ref[i] * v[i];
		}
		return s;
	}

	//operators
	Vector Vector::operator+(void) const
	{
		return *this;
	}
	Vector Vector::operator-(void) const
	{
		return -1 * *this;
	}
	Vector Vector::operator/(double s) const
	{
		return Vector(*this) /= s;
	}

	Vector Vector::operator+(const Vector& v) const
	{
		return Vector(*this) += v;
	}
	Vector Vector::operator-(const Vector& v) const
	{
		return Vector(*this) -= v;
	}

	Vector& Vector::operator=(double s)
	{
		Matrix::operator=(s);
		return *this;
	}
	Vector& Vector::operator=(const double* v)
	{
		Matrix::operator=(v);
		return *this;
	}
	Vector& Vector::operator=(const Vector& v)
	{
		Matrix::operator=(v);
		return *this;
	}
	Vector& Vector::operator=(std::initializer_list<double> list)
	{
		Matrix::operator=(list);
		return *this;
	}

	Vector& Vector::operator+=(double s)
	{
		Matrix::operator+=(s);
		return *this;
	}
	Vector& Vector::operator-=(double s)
	{
		Matrix::operator-=(s);
		return *this;
	}
	Vector& Vector::operator*=(double s)
	{
		Matrix::operator*=(s);
		return *this;
	}
	Vector& Vector::operator/=(double s)
	{
		Matrix::operator/=(s);
		return *this;
	}

	Vector& Vector::operator+=(const double* v)
	{
		Matrix::operator+=(v);
		return *this;
	}
	Vector& Vector::operator-=(const double* v)
	{
		Matrix::operator-=(v);
		return *this;
	}
	Vector& Vector::operator+=(const Vector& v)
	{
		Matrix::operator+=(v);
		return *this;
	}
	Vector& Vector::operator-=(const Vector& v)
	{
		Matrix::operator-=(v);
		return *this;
	}
}