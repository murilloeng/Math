//Math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Solvers/Implicit.hpp"
#include "Math/inc/Solvers/Convergence.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		Convergence::Convergence(Implicit* solver) :
			m_types{1 << uint32_t(Type::Force)}, m_solver{solver}, m_tolerance{1.00e-5}
		{
			return;
		}

		//destructor
		Convergence::~Convergence(void)
		{
			return;
		}

		//data
		double Convergence::tolerance(void) const
		{
			return m_tolerance;
		}
		double Convergence::tolerance(double tolerance)
		{
			return m_tolerance = tolerance;
		}

		uint32_t Convergence::types(void) const
		{
			return m_types;
		}
		uint32_t Convergence::types(uint32_t types)
		{
			return m_types = types;
		}

		uint32_t Convergence::type_set(Type type)
		{
			return m_types = 1 << uint32_t(type);
		}
		uint32_t Convergence::type_add(Type type)
		{
			return m_types |= 1 << uint32_t(type);
		}
		uint32_t Convergence::type_remove(Type type)
		{
			return m_types = m_types & ~(1 << uint32_t(type));
		}

		//check
		bool Convergence::check(void) const
		{
			//data
			const math::Vector r(m_solver->m_r, m_solver->m_size);
			const math::Vector fe(m_solver->m_fe, m_solver->m_size);
			const math::Vector dx(m_solver->m_dx, m_solver->m_size);
			const math::Vector ddx(m_solver->m_ddx, m_solver->m_size);
			//tests
			const bool test_fixed = ~m_types & 1 << uint32_t(Type::Fixed) || r.norm() < m_tolerance;
			const bool test_force = ~m_types & 1 << uint32_t(Type::Force) || r.norm() < m_tolerance * fe.norm();
			const bool test_state = ~m_types & 1 << uint32_t(Type::State) || ddx.norm() < m_tolerance * dx.norm();
			//return
			return test_fixed && test_force && test_state;
		}
	}
}
