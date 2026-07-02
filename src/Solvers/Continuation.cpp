//std
#include <cmath>

//Math
#include "Math/inc/Miscellaneous/util.hpp"
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Solvers/Solver.hpp"
#include "Math/inc/Solvers/Continuation.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		Continuation::Continuation(Solver* solver, Type type) :  m_type{type}, m_solver{solver}
		{
			return;
		}

		//destructor
		Continuation::~Continuation(void)
		{
			return;
		}

		//data
		Continuation::Type Continuation::type(Type type)
		{
			return m_type = type;
		}
		Continuation::Type Continuation::type(void) const
		{
			return m_type;
		}

		//continuation
		double Continuation::predictor(void) const
		{
			//data
			double(Continuation::*pfun[])(void) const = {
				&Continuation::predictor_minimal_norm, 
				&Continuation::predictor_control_load, 
				&Continuation::predictor_control_state, 
				&Continuation::predictor_arc_length_spherical, 
				&Continuation::predictor_arc_length_cylindrical
			};
			//predictor
			for(uint32_t i = 0; 1U << i < uint32_t(Type::Last); i++)
			{
				if(uint32_t(m_type) == 1U << i)
				{
					return (this->*pfun[i])();
				}
			}
			return 0;
		}
		double Continuation::corrector(void) const
		{
			//data
			double(Continuation::*cfun[])(void) const = {
				&Continuation::corrector_minimal_norm, 
				&Continuation::corrector_control_load, 
				&Continuation::corrector_control_state, 
				&Continuation::corrector_arc_length_spherical, 
				&Continuation::corrector_arc_length_cylindrical
			};
			//corrector
			for(uint32_t i = 0; 1U << i < uint32_t(Type::Last); i++)
			{
				if(uint32_t(m_type) == 1U << i)
				{
					return (this->*cfun[i])();
				}
			}
			return 0;
		}

		//types
		double Continuation::predictor_minimal_norm(void) const
		{
			//data
			const math::Vector dx(m_solver->m_dx, m_solver->m_size);
			const math::Vector dxr(m_solver->m_dxr, m_solver->m_size);
			const math::Vector dxt(m_solver->m_dxt, m_solver->m_size);
			//predictor
			return (dx - dxr).inner(dxt) / dxt.inner(dxt);
		}
		double Continuation::corrector_minimal_norm(void) const
		{
			//data
			const math::Vector ddxr(m_solver->m_ddxr, m_solver->m_size);
			const math::Vector ddxt(m_solver->m_ddxt, m_solver->m_size);
			//corrector
			return -ddxr.inner(ddxt) / ddxt.inner(ddxt);
		}
		double Continuation::predictor_control_load(void) const
		{
			return m_solver->m_dp;
		}
		double Continuation::corrector_control_load(void) const
		{
			return 0;
		}
		double Continuation::predictor_control_state(void) const
		{
			const uint32_t index = m_solver->m_watch_dof;
			return (m_solver->m_dx[index] - m_solver->m_dxr[index]) / m_solver->m_dxt[index];
		}
		double Continuation::corrector_control_state(void) const
		{
			const uint32_t index = m_solver->m_watch_dof;
			return -m_solver->m_ddxr[index] / m_solver->m_ddxt[index];
		}
		double Continuation::predictor_arc_length_spherical(void) const
		{
			//data
			const double dl = m_solver->m_dp;
			const math::Vector dx(m_solver->m_dx, m_solver->m_size);
			const math::Vector dxr(m_solver->m_dxr, m_solver->m_size);
			const math::Vector dxt(m_solver->m_dxt, m_solver->m_size);
			//data
			const double b = dxt.inner(dxr);
			const double a = dxt.inner(dxt) + 1;
			const double s = math::sign(dxt.inner(dx));
			const double c = dxr.inner(dxr) - dx.inner(dx) - dl * dl;
			//predictor
			return -b / a + s * sqrt(b * b - c * a) / a;
		}
		double Continuation::corrector_arc_length_spherical(void) const
		{
			//data
			const double dl = m_solver->m_dp;
			const math::Vector dx(m_solver->m_dx, m_solver->m_size);
			const math::Vector ddxr(m_solver->m_ddxr, m_solver->m_size);
			const math::Vector ddxt(m_solver->m_ddxt, m_solver->m_size);
			//data
			const double a = ddxt.inner(ddxt) + 1;
			const double c = ddxr.inner(ddxr + 2 * dx);
			const double b = ddxt.inner(ddxr + dx) + dl;
			const double s = math::sign(ddxt.inner(dx));
			//predictor
			return -b / a + s * sqrt(b * b - c * a) / a;
		}
		double Continuation::predictor_arc_length_cylindrical(void) const
		{
			//data
			const math::Vector dx(m_solver->m_dx, m_solver->m_size);
			const math::Vector dxr(m_solver->m_dxr, m_solver->m_size);
			const math::Vector dxt(m_solver->m_dxt, m_solver->m_size);
			//data
			const double a = dxt.inner(dxt);
			const double b = dxt.inner(dxr);
			const double s = math::sign(dxt.inner(dx));
			const double c = dxr.inner(dxr) - dx.inner(dx);
			//predictor
			return -b / a + s * sqrt(b * b - c * a) / a;
		}
		double Continuation::corrector_arc_length_cylindrical(void) const
		{
			//data
			const math::Vector dx(m_solver->m_dx, m_solver->m_size);
			const math::Vector ddxr(m_solver->m_ddxr, m_solver->m_size);
			const math::Vector ddxt(m_solver->m_ddxt, m_solver->m_size);
			//data
			const double a = ddxt.inner(ddxt);
			const double b = ddxt.inner(ddxr + dx);
			const double c = ddxr.inner(ddxr + 2 * dx);
			const double s = math::sign(ddxt.inner(dx));
			//predictor
			return -b / a + s * sqrt(b * b - c * a) / a;
		}
	}
}