//Math
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/solver.hpp"
#include "Math/Math/inc/solvers/convergence.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		convergence::convergence(void) : 
			m_type(type::force), m_solver(nullptr), m_tolerance(1.00e-5)
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
			const math::vector r(m_solver->m_r, m_solver->m_size);
			const math::vector fe(m_solver->m_fe, m_solver->m_size);
			//return
			return r.norm() < m_tolerance * (m_type == type::fixed ? 1 : fe.norm());
		}
	}
}
