//Math
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/convergence.hpp"
#include "Math/Math/inc/solvers/newton_raphson.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		convergence::convergence(const newton_raphson* solver) : 
			m_type(type::force), m_tolerance(1.00e-5), m_solver(solver)
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
			const math::vector g(m_solver->m_g, m_solver->m_size);
			//return
			return r.norm() < m_tolerance * (m_type == type::fixed ? 1 : g.norm());
		}
	}
}
