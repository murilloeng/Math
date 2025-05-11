//std
#include <cmath>

//Math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/continuation.hpp"

namespace math
{
	namespace solvers
	{
		//constructor
		continuation::continuation(void) : m_type(type::arc_length_cylindrical)
		{
			return;
		}
		continuation::continuation(type type) : m_type(type)
		{
			return;
		}

		//destructor
		continuation::~continuation(void)
		{
			return;
		}

		//continuation
		double continuation::predictor(void) const
		{
			//data
			double(continuation::*pfun[])(void) const = {
				&continuation::predictor_minimal_norm, 
				&continuation::predictor_control_load, 
				&continuation::predictor_control_state, 
				&continuation::predictor_arc_length_spherical, 
				&continuation::predictor_arc_length_cylindrical
			};
			//predictor
			for(uint32_t i = 0; 1U << i < uint32_t(type::last); i++)
			{
				if(uint32_t(m_type) == 1U << i)
				{
					return (this->*pfun[i])();
				}
			}
			return 0;
		}
		double continuation::corrector(void) const
		{
			//data
			double(continuation::*cfun[])(void) const = {
				&continuation::corrector_minimal_norm, 
				&continuation::corrector_control_load, 
				&continuation::corrector_control_state, 
				&continuation::corrector_arc_length_spherical, 
				&continuation::corrector_arc_length_cylindrical
			};
			//corrector
			for(uint32_t i = 0; 1U << i < uint32_t(type::last); i++)
			{
				if(uint32_t(m_type) == 1U << i)
				{
					return (this->*cfun[i])();
				}
			}
			return 0;
		}

		//types
		double continuation::predictor_minimal_norm(void) const
		{
			//data
			const math::vector dx(m_dx, m_size);
			const math::vector dxr(m_dxr, m_size);
			const math::vector dxt(m_dxt, m_size);
			//predictor
			return (dx - dxr).inner(dxt) / dxt.inner(dxt);
		}
		double continuation::corrector_minimal_norm(void) const
		{
			//data
			const math::vector ddxr(m_ddxr, m_size);
			const math::vector ddxt(m_ddxt, m_size);
			//corrector
			return -ddxr.inner(ddxt) / ddxt.inner(ddxt);
		}
		double continuation::predictor_control_load(void) const
		{
			return *m_dp;
		}
		double continuation::corrector_control_load(void) const
		{
			return 0;
		}
		double continuation::predictor_control_state(void) const
		{
			return (m_dx[m_index] - m_dxr[m_index]) / m_dxt[m_index];
		}
		double continuation::corrector_control_state(void) const
		{
			return -m_ddxr[m_index] / m_ddxt[m_index];
		}
		double continuation::predictor_arc_length_spherical(void) const
		{
			//data
			const double dl = *m_dp;
			const math::vector dx(m_dx, m_size);
			const math::vector dxr(m_dxr, m_size);
			const math::vector dxt(m_dxt, m_size);
			//data
			const double b = dxt.inner(dxr);
			const double a = dxt.inner(dxt) + 1;
			const double s = math::sign(dxt.inner(dx));
			const double c = dxr.inner(dxr) - dx.inner(dx) - dl * dl;
			//predictor
			return -b / a + s * sqrt(b * b - c * a) / a;
		}
		double continuation::corrector_arc_length_spherical(void) const
		{
			//data
			const double dl = *m_dp;
			const math::vector dx(m_dx, m_size);
			const math::vector ddxr(m_ddxr, m_size);
			const math::vector ddxt(m_ddxt, m_size);
			//data
			const double a = ddxt.inner(ddxt) + 1;
			const double c = ddxr.inner(ddxr + 2 * dx);
			const double b = ddxt.inner(ddxr + dx) + dl;
			const double s = math::sign(ddxt.inner(dx));
			//predictor
			return -b / a + s * sqrt(b * b - c * a) / a;
		}
		double continuation::predictor_arc_length_cylindrical(void) const
		{
			//data
			const math::vector dx(m_dx, m_size);
			const math::vector dxr(m_dxr, m_size);
			const math::vector dxt(m_dxt, m_size);
			//data
			const double a = dxt.inner(dxt);
			const double b = dxt.inner(dxr);
			const double s = math::sign(dxt.inner(dx));
			const double c = dxr.inner(dxr) - dx.inner(dx);
			//predictor
			return -b / a + s * sqrt(b * b - c * a) / a;
		}
		double continuation::corrector_arc_length_cylindrical(void) const
		{
			//data
			const math::vector dx(m_dx, m_size);
			const math::vector ddxr(m_ddxr, m_size);
			const math::vector ddxt(m_ddxt, m_size);
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