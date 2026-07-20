//Math
#include "Math/inc/Eigen/Eigen.hpp"

namespace math
{
	namespace eigen
	{
		//constructor
		Eigen::Eigen(uint32_t order) : m_order{order}, m_tolerance{0}
		{
			return;
		}
		
		//destructor
		Eigen::~Eigen(void)
		{
			return;
		}

		//data
		uint32_t Eigen::order(void) const
		{
			return m_order;
		}
		uint32_t Eigen::order(uint32_t order)
		{
			return m_order = order;
		}

		double Eigen::tolerance(void) const
		{
			return m_tolerance;
		}
		double Eigen::tolerance(double tolerance)
		{
			return m_tolerance = tolerance;
		}
	}
}