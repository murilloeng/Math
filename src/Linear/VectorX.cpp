//std
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Linear/VectorX.hpp"
#include "Math/inc/Linear/MatrixX.hpp"

namespace math
{
	//constructors
	VectorX::VectorX(uint32_t rows) : m_own{true}, m_rows{rows}
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
	VectorX::VectorX(const VectorX& v) : VectorX(v.m_rows)
	{
		memcpy(m_data_ptr, v.m_data_ref, m_rows * sizeof(double));
	}
	VectorX::VectorX(double* ptr, uint32_t rows) : m_own{false}, m_rows{rows}, m_data_ptr{ptr}, m_data_ref{ptr}
	{
		return;
	}
	VectorX::VectorX(const double* ptr, uint32_t rows) : m_own{false}, m_rows{rows}, m_data_ptr{nullptr}, m_data_ref{ptr}
	{
		return;
	}
	VectorX::VectorX(std::initializer_list<double> list) : VectorX((uint32_t) list.size())
	{
		memcpy(m_data_ptr, list.begin(), list.size() * sizeof(double));
	}

	//destructor
	VectorX::~VectorX(void)
	{
		if(m_own && m_data_ptr != m_data_mem)
		{
			delete[] m_data_ptr;
		}
	}

	//size
	uint32_t VectorX::rows(void) const
	{
		return m_rows;
	}

	//data
	double* VectorX::data(void)
	{
		return m_data_ptr;
	}
	const double* VectorX::data(void) const
	{
		return m_data_ref;
	}

	//direction
	double VectorX::norm(void) const
	{
		double s = 0;
		for(uint32_t i = 0; i < m_rows; i++)
		{
			s += m_data_ref[i] * m_data_ref[i];
		}
		return sqrt(s);
	}
	VectorX& VectorX::normalize(void)
	{
		return *this /= norm();
	}
	VectorX VectorX::unit(void) const
	{
		return *this / norm();
	}

	//utility
	VectorX& VectorX::zeros(void)
	{
		memset(m_data_ptr, 0, m_rows * sizeof(double));
		return *this;
	}
	VectorX& VectorX::randu(double a, double b)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] = a + (b - a) * rand() / RAND_MAX;
		}
		return *this;
	}
	void VectorX::print(const char* header, double cutoff) const
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
	MatrixX VectorX::outer(void) const
	{
		return outer(*this);
	}
	MatrixX VectorX::outer(const VectorX& v) const
	{
		MatrixX m(m_rows, v.m_rows);
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < v.m_rows; j++)
			{
				m(i, j) = m_data_ref[i] * v.m_data_ref[j];
			}
		}
		return m;
	}
	double VectorX::inner(const VectorX& v) const
	{
		double s = 0;
		if(m_rows != v.m_rows)
		{
			throw std::runtime_error("VectorX inner product called with inconsistent dimensions!");
		}
		for(uint32_t i = 0; i < m_rows; i++)
		{
			s += m_data_ref[i] * v.m_data_ref[i];
		}
		return s;
	}
	double VectorX::inner(const double* v) const
	{
		double s = 0;
		for(uint32_t i = 0; i < m_rows; i++)
		{
			s += m_data_ref[i] * v[i];
		}
		return s;
	}

	//operators
	VectorX VectorX::operator+(void) const
	{
		return *this;
	}
	VectorX VectorX::operator-(void) const
	{
		return -1 * VectorX(*this);
	}
	VectorX VectorX::operator/(double s) const
	{
		return VectorX(*this) /= s;
	}

	VectorX VectorX::operator+(const VectorX& v) const
	{
		return VectorX(*this) += v;
	}
	VectorX VectorX::operator-(const VectorX& v) const
	{
		return VectorX(*this) -= v;
	}

	VectorX& VectorX::operator=(double s)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] = s;
		}
		return *this;
	}
	VectorX& VectorX::operator=(const double* v)
	{
		memcpy(m_data_ptr, v, m_rows * sizeof(double));
		return *this;
	}
	VectorX& VectorX::operator=(const VectorX& v)
	{
		//check
		if(m_rows != v.m_rows)
		{
			throw std::runtime_error("VectorX assign has incompatible dimensions!");
		}
		//assign
		memcpy(m_data_ptr, v.m_data_ref, m_rows * sizeof(double));
		//return
		return *this;
	}
	VectorX& VectorX::operator=(std::initializer_list<double> list)
	{
		//check
		if(list.size() != m_rows)
		{
			throw std::runtime_error("VectorX assign has incompatible dimensions!");
		}
		//assign
		memcpy(m_data_ptr, list.begin(), m_rows * sizeof(double));
		//return
		return *this;
	}

	VectorX& VectorX::operator+=(double s)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] += s;
		}
		return *this;
	}
	VectorX& VectorX::operator-=(double s)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] -= s;
		}
		return *this;
	}
	VectorX& VectorX::operator*=(double s)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] *= s;
		}
		return *this;
	}
	VectorX& VectorX::operator/=(double s)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] /= s;
		}
		return *this;
	}

	VectorX& VectorX::operator+=(const double* v)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] += v[i];
		}
		return *this;
	}
	VectorX& VectorX::operator-=(const double* v)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] -= v[i];
		}
		return *this;
	}
	VectorX& VectorX::operator+=(const VectorX& v)
	{
		if(m_rows != v.m_rows)
		{
			throw std::runtime_error("VectorX increment has incompatible dimensions!");
		}
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] += v.m_data_ref[i];
		}
		return *this;
	}
	VectorX& VectorX::operator-=(const VectorX& v)
	{
		if(m_rows != v.m_rows)
		{
			throw std::runtime_error("VectorX increment has incompatible dimensions!");
		}
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] -= v.m_data_ref[i];
		}
		return *this;
	}

	double& VectorX::operator[](uint32_t index)
	{
		return m_data_ptr[index];
	}
	const double& VectorX::operator[](uint32_t index) const
	{
		return m_data_ref[index];
	}

	VectorX operator*(double s, const VectorX& v)
	{
		return VectorX(v) *= s;
	}
}