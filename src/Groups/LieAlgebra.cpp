//std
#include <stdexcept>

//Math
#include "Math/inc/Groups/LieAlgebra.hpp"

namespace math
{
	//constructor
	LieAlgebra::LieAlgebra(void) : m_dimension{0}
	{
		return;
	}

	//destructor
	LieAlgebra::~LieAlgebra(void)
	{
		return;
	}

	//data
	double* LieAlgebra::data(void)
	{
		return m_data;
	}
	const double* LieAlgebra::data(void) const
	{
		return m_data;
	}
	uint32_t LieAlgebra::dimension(void) const
	{
		return m_dimension;
	}

	//operators
	LieAlgebra& LieAlgebra::operator*=(double s)
	{
		for(uint32_t i = 0; i < m_dimension; i++) m_data[i] *= s;
		return *this;
	}
	LieAlgebra& LieAlgebra::operator/=(double s)
	{
		for(uint32_t i = 0; i < m_dimension; i++) m_data[i] /= s;
		return *this;
	}
	LieAlgebra& LieAlgebra::operator+=(const LieAlgebra& lie_algebra)
	{
		//check
		if(lie_algebra.m_dimension != m_dimension)
		{
			throw std::runtime_error("Lie algebra operation with incompatible dimensions!");
		}
		//compute
		for(uint32_t i = 0; i < m_dimension; i++) m_data[i] += lie_algebra.m_data[i];
		//return
		return *this;
	}
	LieAlgebra& LieAlgebra::operator-=(const LieAlgebra& lie_algebra)
	{
		//check
		if(lie_algebra.m_dimension != m_dimension)
		{
			throw std::runtime_error("Lie algebra operation with incompatible dimensions!");
		}
		//compute
		for(uint32_t i = 0; i < m_dimension; i++) m_data[i] -= lie_algebra.m_data[i];
		//return
		return *this;
	}

	LieAlgebra LieAlgebra::operator+(void) const
	{
		return *this;
	}
	LieAlgebra LieAlgebra::operator-(void) const
	{
		return LieAlgebra(*this) *= -1;
	}
	LieAlgebra LieAlgebra::operator*(double s) const
	{
		return LieAlgebra(*this) *= s;
	}
	LieAlgebra LieAlgebra::operator/(double s) const
	{
		return LieAlgebra(*this) /= s;
	}
	LieAlgebra LieAlgebra::operator+(const LieAlgebra& lie_algebra) const
	{
		return LieAlgebra(*this) += lie_algebra;
	}
	LieAlgebra LieAlgebra::operator-(const LieAlgebra& lie_algebra) const
	{
		return LieAlgebra(*this) -= lie_algebra;
	}

	LieAlgebra operator*(double s, const LieAlgebra& lie_algebra)
	{
		return LieAlgebra(lie_algebra) *= s;
	}
}