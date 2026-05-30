//std
#include <cmath>

//Math
#include "Math/inc/misc/util.hpp"
#include "Math/inc/linear/vector.hpp"
#include "Math/inc/solvers/Solver.hpp"
#include "Math/inc/solvers/Continuation.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		Continuation::Continuation(void) : m_type(Type::ArcLengthCylindrical)
		{
			return;
		}
		Continuation::Continuation(Type type) : m_type(type)
		{
			return;
		}

		//destructor
		Continuation::~Continuation(void)
		{
			return;
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
			const math::vector dx(m_solver->m_dx, m_solver->m_size);
			const math::vector dxr(m_solver->m_dxr, m_solver->m_size);
			const math::vector dxt(m_solver->m_dxt, m_solver->m_size);
			//predictor
			return (dx - dxr).inner(dxt) / dxt.inner(dxt);
		}
		double Continuation::corrector_minimal_norm(void) const
		{
			//data
			const math::vector ddxr(m_solver->m_ddxr, m_solver->m_size);
			const math::vector ddxt(m_solver->m_ddxt, m_solver->m_size);
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
			const math::vector dx(m_solver->m_dx, m_solver->m_size);
			const math::vector dxr(m_solver->m_dxr, m_solver->m_size);
			const math::vector dxt(m_solver->m_dxt, m_solver->m_size);
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
			const math::vector dx(m_solver->m_dx, m_solver->m_size);
			const math::vector ddxr(m_solver->m_ddxr, m_solver->m_size);
			const math::vector ddxt(m_solver->m_ddxt, m_solver->m_size);
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
			const math::vector dx(m_solver->m_dx, m_solver->m_size);
			const math::vector dxr(m_solver->m_dxr, m_solver->m_size);
			const math::vector dxt(m_solver->m_dxt, m_solver->m_size);
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
			const math::vector dx(m_solver->m_dx, m_solver->m_size);
			const math::vector ddxr(m_solver->m_ddxr, m_solver->m_size);
			const math::vector ddxt(m_solver->m_ddxt, m_solver->m_size);
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