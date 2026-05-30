//Math
#include "Math/inc/linear/vector.hpp"
#include "Math/inc/solvers/Solver.hpp"
#include "Math/inc/solvers/Convergence.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		Convergence::Convergence(void) : 
			m_type(Type::Force), m_solver(nullptr), m_tolerance(1.00e-5)
		{
			return;
		}
		
		//destructor
		Convergence::~Convergence(void)
		{
			return;
		}

		//check
		bool Convergence::check(void) const
		{
			//data
			const math::vector r(m_solver->m_r, m_solver->m_size);
			const math::vector fe(m_solver->m_fe, m_solver->m_size);
			//return
			return r.norm() < m_tolerance * (m_type == Type::Fixed ? 1 : fe.norm());
		}
	}
}
