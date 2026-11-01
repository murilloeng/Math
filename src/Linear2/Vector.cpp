//std
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Linear2/Vector.hpp"

namespace math
{
	namespace linear
	{
		//constructors
		Vector::Vector(uint32_t rows) : Matrix{rows, 1}
		{
			return;
		}
		Vector::Vector(const Vector& v) : Matrix{v}
		{
			return;
		}
		Vector::Vector(const Matrix& m) : Matrix{m}
		{
			if(m_cols != 1)
			{
				throw std::runtime_error("Vector constructor has incompatible dimensions!");
			}
		}
		Vector::Vector(double* ptr, uint32_t rows) : Matrix(ptr, rows, 1)
		{
			return;
		}
		Vector::Vector(const double* ptr, uint32_t rows) : Matrix{ptr, rows, 1}
		{
			return;
		}
		Vector::Vector(std::initializer_list<double> list) : Matrix{(uint32_t) list.size(), 1}
		{
			memcpy(m_data_ptr, std::data(list), list.size() * sizeof(double));
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
			if(m_rows != v.m_rows)
			{
				throw std::runtime_error("Vector inner product called with inconsistent dimensions!");
			}
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
}