//std
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>

//math
#include "Math/inc/linear/vector.hpp"

namespace math
{
	//constructors
	vector::vector(void) : matrix(0, 1)
	{
		return;
	}
	vector::vector(uint32_t rows, mode init) : matrix(rows, 1, init)
	{
		return;
	}
	vector::vector(const matrix& m) : matrix(m.rows() * m.cols(), 1)
	{
		memcpy(m_ptr, m.data(), m_rows * sizeof(double));
	}
	vector::vector(const double* ptr, uint32_t rows) : matrix(ptr, rows, 1)
	{
		return;
	}
	vector::vector(std::initializer_list<double> list) : matrix(list.size(), 1)
	{
		memcpy(m_ptr, std::data(list), list.size() * sizeof(double));
	}
	vector::vector(double* ptr, uint32_t rows, mode init) : matrix(ptr, rows, 1, init)
	{
		return;
	}

	//destructor
	vector::~vector(void)
	{
		return;
	}

	//size
	vector& vector::resize(uint32_t rows)
	{
		matrix::resize(rows, 1);
		return *this;
	}
	vector& vector::resize(uint32_t rows, uint32_t cols)
	{
		//check
		if(cols != 1)
		{
			fprintf(stderr, "Error: Resize called on vector with number of columns not equal to 1!");
			exit(EXIT_FAILURE);
		}
		//resize
		matrix::resize(rows, 1);
		return *this;
	}

	//linear
	vector vector::unit(void) const
	{
		return *this / norm();
	}
	matrix vector::outer(void) const
	{
		return outer(*this);
	}
	matrix vector::outer(const vector& v) const
	{
		matrix m(m_rows, v.m_rows);
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < v.m_rows; j++)
			{
				m(i, j) = m_ref[i] * v.m_ref[j];
			}
		}
		return m;
	}
	double vector::inner(const vector& v) const
	{
		double s = 0;
		for(uint32_t i = 0; i < m_rows; i++)
		{
			s += m_ref[i] * v.m_ref[i];
		}
		return s;
	}
	double vector::inner(const double* v) const
	{
		double s = 0;
		for(uint32_t i = 0; i < m_rows; i++)
		{
			s += m_ref[i] * v[i];
		}
		return s;
	}

	//operators
	vector vector::operator+(void) const
	{
		return *this;
	}
	vector vector::operator-(void) const
	{
		return -1 * *this;
	}
	vector vector::operator/(double s) const
	{
		return vector(*this) /= s;
	}

	vector vector::operator+(const vector& v) const
	{
		return vector(*this) += v;
	}
	vector vector::operator-(const vector& v) const
	{
		return vector(*this) -= v;
	}

	vector& vector::operator=(double s)
	{
		matrix::operator=(s);
		return *this;
	}
	vector& vector::operator=(const double* v)
	{
		matrix::operator=(v);
		return *this;
	}
	vector& vector::operator=(const vector& v)
	{
		matrix::operator=(v);
		return *this;
	}
	vector& vector::operator=(std::initializer_list<double> list)
	{
		matrix::operator=(list);
		return *this;
	}

	vector& vector::operator+=(double s)
	{
		matrix::operator+=(s);
		return *this;
	}
	vector& vector::operator-=(double s)
	{
		matrix::operator-=(s);
		return *this;
	}
	vector& vector::operator*=(double s)
	{
		matrix::operator*=(s);
		return *this;
	}
	vector& vector::operator/=(double s)
	{
		matrix::operator/=(s);
		return *this;
	}

	vector& vector::operator+=(const double* v)
	{
		matrix::operator+=(v);
		return *this;
	}
	vector& vector::operator-=(const double* v)
	{
		matrix::operator-=(v);
		return *this;
	}
	vector& vector::operator+=(const vector& v)
	{
		matrix::operator+=(v);
		return *this;
	}
	vector& vector::operator-=(const vector& v)
	{
		matrix::operator-=(v);
		return *this;
	}
}