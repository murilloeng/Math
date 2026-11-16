//std
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Linear/VectorN.hpp"
#include "Math/inc/Linear/MatrixN.hpp"

namespace math
{
	//constructors
	VectorN::VectorN(uint32_t rows) : m_own{true}, m_rows{rows}
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
	VectorN::VectorN(const VectorN& v) : VectorN(v.m_rows)
	{
		memcpy(m_data_ptr, v.m_data_ref, m_rows * sizeof(double));
	}
	VectorN::VectorN(double* ptr, uint32_t rows) : m_own{false}, m_rows{rows}, m_data_ptr{ptr}, m_data_ref{ptr}
	{
		return;
	}
	VectorN::VectorN(const double* ptr, uint32_t rows) : m_own{false}, m_rows{rows}, m_data_ptr{nullptr}, m_data_ref{ptr}
	{
		return;
	}
	VectorN::VectorN(std::initializer_list<double> list) : VectorN((uint32_t) list.size())
	{
		memcpy(m_data_ptr, list.begin(), list.size() * sizeof(double));
	}

	//destructor
	VectorN::~VectorN(void)
	{
		if(m_own && m_data_ptr != m_data_mem)
		{
			delete[] m_data_ptr;
		}
	}

	//size
	uint32_t VectorN::rows(void) const
	{
		return m_rows;
	}

	//data
	double* VectorN::data(void)
	{
		return m_data_ptr;
	}
	const double* VectorN::data(void) const
	{
		return m_data_ref;
	}

	//direction
	double VectorN::norm(void) const
	{
		double s = 0;
		for(uint32_t i = 0; i < m_rows; i++)
		{
			s += m_data_ref[i] * m_data_ref[i];
		}
		return sqrt(s);
	}
	VectorN& VectorN::normalize(void)
	{
		return *this /= norm();
	}
	VectorN VectorN::unit(void) const
	{
		return *this / norm();
	}

	//utility
	VectorN& VectorN::zeros(void)
	{
		memset(m_data_ptr, 0, m_rows * sizeof(double));
		return *this;
	}
	VectorN& VectorN::randu(double a, double b)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] = a + (b - a) * rand() / RAND_MAX;
		}
		return *this;
	}
	void VectorN::print(const char* header, double cutoff) const
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
	MatrixN VectorN::outer(void) const
	{
		return outer(*this);
	}
	MatrixN VectorN::outer(const VectorN& v) const
	{
		MatrixN m(m_rows, v.m_rows);
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < v.m_rows; j++)
			{
				m(i, j) = m_data_ref[i] * v.m_data_ref[j];
			}
		}
		return m;
	}
	double VectorN::inner(const VectorN& v) const
	{
		double s = 0;
		if(m_rows != v.m_rows)
		{
			throw std::runtime_error("VectorN inner product called with inconsistent dimensions!");
		}
		for(uint32_t i = 0; i < m_rows; i++)
		{
			s += m_data_ref[i] * v.m_data_ref[i];
		}
		return s;
	}
	double VectorN::inner(const double* v) const
	{
		double s = 0;
		for(uint32_t i = 0; i < m_rows; i++)
		{
			s += m_data_ref[i] * v[i];
		}
		return s;
	}

	//operators
	VectorN VectorN::operator+(void) const
	{
		return *this;
	}
	VectorN VectorN::operator-(void) const
	{
		return -1 * VectorN(*this);
	}
	VectorN VectorN::operator/(double s) const
	{
		return VectorN(*this) /= s;
	}

	VectorN VectorN::operator+(const VectorN& v) const
	{
		return VectorN(*this) += v;
	}
	VectorN VectorN::operator-(const VectorN& v) const
	{
		return VectorN(*this) -= v;
	}

	VectorN& VectorN::operator=(double s)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] = s;
		}
		return *this;
	}
	VectorN& VectorN::operator=(const double* v)
	{
		memcpy(m_data_ptr, v, m_rows * sizeof(double));
		return *this;
	}
	VectorN& VectorN::operator=(const VectorN& v)
	{
		//check
		if(m_rows != v.m_rows)
		{
			throw std::runtime_error("VectorN assign has incompatible dimensions!");
		}
		//assign
		memcpy(m_data_ptr, v.m_data_ref, m_rows * sizeof(double));
		//return
		return *this;
	}
	VectorN& VectorN::operator=(std::initializer_list<double> list)
	{
		//check
		if(list.size() != m_rows)
		{
			throw std::runtime_error("VectorN assign has incompatible dimensions!");
		}
		//assign
		memcpy(m_data_ptr, list.begin(), m_rows * sizeof(double));
		//return
		return *this;
	}

	VectorN& VectorN::operator+=(double s)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] += s;
		}
		return *this;
	}
	VectorN& VectorN::operator-=(double s)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] -= s;
		}
		return *this;
	}
	VectorN& VectorN::operator*=(double s)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] *= s;
		}
		return *this;
	}
	VectorN& VectorN::operator/=(double s)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] /= s;
		}
		return *this;
	}

	VectorN& VectorN::operator+=(const double* v)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] += v[i];
		}
		return *this;
	}
	VectorN& VectorN::operator-=(const double* v)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] -= v[i];
		}
		return *this;
	}
	VectorN& VectorN::operator+=(const VectorN& v)
	{
		if(m_rows != v.m_rows)
		{
			throw std::runtime_error("VectorN increment has incompatible dimensions!");
		}
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] += v.m_data_ref[i];
		}
		return *this;
	}
	VectorN& VectorN::operator-=(const VectorN& v)
	{
		if(m_rows != v.m_rows)
		{
			throw std::runtime_error("VectorN increment has incompatible dimensions!");
		}
		for(uint32_t i = 0; i < m_rows; i++)
		{
			m_data_ptr[i] -= v.m_data_ref[i];
		}
		return *this;
	}

	double& VectorN::operator[](uint32_t index)
	{
		return m_data_ptr[index];
	}
	const double& VectorN::operator[](uint32_t index) const
	{
		return m_data_ref[index];
	}

	VectorN operator*(double s, const VectorN& v)
	{
		return VectorN(v) *= s;
	}
}