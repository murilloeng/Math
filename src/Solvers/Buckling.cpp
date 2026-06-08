//Math
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
			m_stiffness_0(m_K);
			math::Vector F(m_fe, m_size);
			math::Vector x(m_x_new, m_size);
			math::Matrix K(m_K, m_size, m_size);
			//linear
			K.solve(x, F);
			m_stiffness_d(m_M, m_x_new);
			
		}
	}
}