//std
#include <cstdio>

//Math
#include "Math/inc/Linear/Eigen.hpp"
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Solvers/Buckling.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		Buckling::Buckling(void)
		{
			return;
		}
		
		//destructor
		Buckling::~Buckling(void)
		{
			return;
		}

		//data
		uint32_t Buckling::state_set(void) const
		{
			return uint32_t(State::x);
		}
		uint32_t Buckling::force_set(void) const
		{
			return 0;
		}
		uint32_t Buckling::tangent_set(void) const
		{
			return uint32_t(Tangent::K) | uint32_t(Tangent::M);
		}

		//solve
		void Buckling::solve(void)
		{
			//data
			Eigen eigen;
			eigen.data(0, m_M);
			eigen.data(1, m_K);
			eigen.index_min(1);
			eigen.index_max(1);
			eigen.type(Eigen::Type::Index);
			//linear
			m_stiffness(m_K, m_x_new);
			Matrix(m_K, m_size, m_size).solve(m_dx, m_fe);
			//apply
			apply();
			m_stiffness(m_M, m_x_new);
			for(uint32_t i = 0; i < m_size * m_size; i++)
			{
				m_M[i] = m_K[i] - m_M[i];
			}
			//compute
			eigen.compute();
			printf("Load: %+.6e\n", eigen.eigenvalue(0, 0));
		}
	}
}