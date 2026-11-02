//std
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Linear2/Vector.hpp"
#include "Math/inc/Linear2/Matrix.hpp"

namespace math
{
	namespace linear
	{
		//constructors
		Vector::Vector(uint32_t rows) : m_own{true}, m_rows{rows}
		{
			if(m_rows <= VECTOR_STATIC_SIZE)
			{
				m_data_ref = m_data_ptr = m_data_mem;
			}
			else
			{
				m_data_ref = m_data_ptr = new double[m_rows];
			}
		}
		Vector::Vector(const Vector& v) : Vector(v.m_rows)
		{
			memcpy(m_data_ptr, v.m_data_ref, m_rows * sizeof(double));
		}
		Vector::Vector(double* ptr, uint32_t rows) : m_own{false}, m_rows{rows}, m_data_ptr{ptr}, m_data_ref{ptr}
		{
			return;
		}
		Vector::Vector(const double* ptr, uint32_t rows) : m_own{false}, m_rows{rows}, m_data_ptr{nullptr}, m_data_ref{ptr}
		{
			return;
		}
		Vector::Vector(std::initializer_list<double> list) : Vector((uint32_t) list.size())
		{
			memcpy(m_data_ptr, list.begin(), list.size() * sizeof(double));
		}

		//destructor
		Vector::~Vector(void)
		{
			if(m_own && m_data_ptr != m_data_mem)
			{
				delete[] m_data_ptr;
			}
		}

		//size
		uint32_t Vector::rows(void) const
		{
			return m_rows;
		}

		//data
		double* Vector::data(void)
		{
			return m_data_ptr;
		}
		const double* Vector::data(void) const
		{
			return m_data_ref;
		}

		//direction
		double Vector::norm(void) const
		{
			double s = 0;
			for(uint32_t i = 0; i < m_rows; i++)
			{
				s += m_data_ref[i] * m_data_ref[i];
			}
			return sqrt(s);
		}
		Vector& Vector::normalize(void)
		{
			return *this /= norm();
		}
		Vector Vector::unit(void) const
		{
			return *this / norm();
		}

		//utility
		Vector& Vector::zeros(void)
		{
			memset(m_data_ptr, 0, m_rows * sizeof(double));
			return *this;
		}
		Vector& Vector::randu(double a, double b)
		{
			for(uint32_t i = 0; i < m_rows; i++)
			{
				m_data_ptr[i] = a + (b - a) * rand() / RAND_MAX;
			}
			return *this;
		}
		void Vector::print(const char* header, double cutoff) const
		{
			if(strlen(header) != 0)
			{
				printf("%s\n", header);
			}
			for(uint32_t i = 0; i < m_rows; i++)
			{
				if(cutoff != 0 && fabs(m_data_ref[i]) < cutoff)
				{
					printf("---------\n");
				}
				else
				{
					printf("%+.2e\n", m_data_ref[i]);
				}
			}
		}

		//products
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
			return -1 * Vector(*this);
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
			for(uint32_t i = 0; i < m_rows; i++)
			{
				m_data_ptr[i] = s;
			}
			return *this;
		}
		Vector& Vector::operator=(const double* v)
		{
			memcpy(m_data_ptr, v, m_rows * sizeof(double));
			return *this;
		}
		Vector& Vector::operator=(const Vector& v)
		{
			//check
			if(m_rows != v.m_rows)
			{
				throw std::runtime_error("Vector assign has incompatible dimensions!");
			}
			//assign
			memcpy(m_data_ptr, v.m_data_ref, m_rows * sizeof(double));
			//return
			return *this;
		}
		Vector& Vector::operator=(std::initializer_list<double> list)
		{
			//check
			if(list.size() != m_rows)
			{
				throw std::runtime_error("Vector assign has incompatible dimensions!");
			}
			//assign
			memcpy(m_data_ptr, list.begin(), m_rows * sizeof(double));
			//return
			return *this;
		}

		Vector& Vector::operator+=(double s)
		{
			for(uint32_t i = 0; i < m_rows; i++)
			{
				m_data_ptr[i] += s;
			}
			return *this;
		}
		Vector& Vector::operator-=(double s)
		{
			for(uint32_t i = 0; i < m_rows; i++)
			{
				m_data_ptr[i] -= s;
			}
			return *this;
		}
		Vector& Vector::operator*=(double s)
		{
			for(uint32_t i = 0; i < m_rows; i++)
			{
				m_data_ptr[i] *= s;
			}
			return *this;
		}
		Vector& Vector::operator/=(double s)
		{
			for(uint32_t i = 0; i < m_rows; i++)
			{
				m_data_ptr[i] /= s;
			}
			return *this;
		}

		Vector& Vector::operator+=(const double* v)
		{
			for(uint32_t i = 0; i < m_rows; i++)
			{
				m_data_ptr[i] += v[i];
			}
			return *this;
		}
		Vector& Vector::operator-=(const double* v)
		{
			for(uint32_t i = 0; i < m_rows; i++)
			{
				m_data_ptr[i] -= v[i];
			}
			return *this;
		}
		Vector& Vector::operator+=(const Vector& v)
		{
			if(m_rows != v.m_rows)
			{
				throw std::runtime_error("Vector increment has incompatible dimensions!");
			}
			for(uint32_t i = 0; i < m_rows; i++)
			{
				m_data_ptr[i] += v.m_data_ref[i];
			}
			return *this;
		}
		Vector& Vector::operator-=(const Vector& v)
		{
			if(m_rows != v.m_rows)
			{
				throw std::runtime_error("Vector increment has incompatible dimensions!");
			}
			for(uint32_t i = 0; i < m_rows; i++)
			{
				m_data_ptr[i] -= v.m_data_ref[i];
			}
			return *this;
		}

		double& Vector::operator[](uint32_t index)
		{
			return m_data_ptr[index];
		}
		const double& Vector::operator[](uint32_t index) const
		{
			return m_data_ref[index];
		}

		Vector operator*(double s, const Vector& v)
		{
			return Vector(v) *= s;
		}
	}
}