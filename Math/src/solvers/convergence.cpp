//Math
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/convergence.hpp"
#include "Math/Math/inc/solvers/newton_raphson.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		convergence::convergence(void) : 
			m_type(type::force), m_r(nullptr), m_g(nullptr), m_size(0), m_tolerance(1.00e-5)
		{
			return;
		}
		
		//destructor
		convergence::~convergence(void)
		{
			return;
		}

		//check
		bool convergence::check(void) const
		{
			//data
			const math::vector r(m_r, m_size);
			const math::vector g(m_g, m_size);
			//return
			return r.norm() < m_tolerance * (m_type == type::fixed ? 1 : g.norm());
		}
	}
}
