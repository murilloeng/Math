//std
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>

//math
#include "inc/linear/vector.hpp"

namespace math
{
	//constructors
	vector::vector(void) : matrix(0, 1)
	{
		return;
	}
	vector::vector(unsigned rows, mode init) : matrix(rows, 1, init)
	{
		return;
	}
	vector::vector(const matrix& m) : matrix(m.rows() * m.cols(), 1)
	{
		memcpy(m_ptr, m.data(), m_rows * sizeof(double));
	}
	vector::vector(const double* ptr, unsigned rows) : matrix(ptr, rows, 1)
	{
		return;
	}
	vector::vector(std::initializer_list<double> list) : matrix(list.size(), 1)
	{
		memcpy(m_ptr, std::data(list), list.size() * sizeof(double));
	}
	vector::vector(double* ptr, unsigned rows, mode init) : matrix(ptr, rows, 1, init)
	{
		return;
	}

	//destructor
	vector::~vector(void)
	{
		return;
	}

	//size
	vector& vector::resize(unsigned rows)
	{
		matrix::resize(rows, 1);
		return *this;
	}
	vector& vector::resize(unsigned rows, unsigned cols)
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
		for(unsigned i = 0; i < m_rows; i++)
		{
			for(unsigned j = 0; j < v.m_rows; j++)
			{
				m(i, j) = m_ref[i] * v.m_ref[j];
			}
		}
		return m;
	}
	double vector::inner(const vector& v) const
	{
		double s = 0;
		for(unsigned i = 0; i < m_rows; i++)
		{
			s += m_ref[i] * v.m_ref[i];
		}
		return s;
	}
}