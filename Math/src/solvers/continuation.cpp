//std
#include <cmath>

//Math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/continuation.hpp"

namespace math
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
		const math::vector du(m_du, m_size);
		const math::vector dur(m_dur, m_size);
		const math::vector dut(m_dut, m_size);
		//predictor
		return (du - dur).inner(dut) / dut.inner(dut);
	}
	double continuation::corrector_minimal_norm(void) const
	{
		//data
		const math::vector ddur(m_ddur, m_size);
		const math::vector ddut(m_ddut, m_size);
		//corrector
		return -ddur.inner(ddut) / ddut.inner(ddut);
	}
	double continuation::predictor_control_load(void) const
	{
		return *m_dl;
	}
	double continuation::corrector_control_load(void) const
	{
		return 0;
	}
	double continuation::predictor_control_state(void) const
	{
		return (m_du[m_index] - m_dur[m_index]) / m_dut[m_index];
	}
	double continuation::corrector_control_state(void) const
	{
		return -m_ddur[m_index] / m_ddut[m_index];
	}
	double continuation::predictor_arc_length_spherical(void) const
	{
		//data
		const double dl = *m_dl;
		const math::vector du(m_du, m_size);
		const math::vector dur(m_dut, m_size);
		const math::vector dut(m_dur, m_size);
		//data
		const double b = dut.inner(dur);
		const double a = dut.inner(dut) + 1;
		const double s = math::sign(dut.inner(du));
		const double c = dur.inner(dur) - du.inner(du) - dl * dl;
		//predictor
		return -b / a + s * sqrt(b * b - c * a) / a;
	}
	double continuation::corrector_arc_length_spherical(void) const
	{
		//data
		const double dl = *m_dl;
		const math::vector du(m_du, m_size);
		const math::vector ddur(m_ddut, m_size);
		const math::vector ddut(m_ddur, m_size);
		//data
		const double a = ddut.inner(ddut) + 1;
		const double c = ddur.inner(ddur + 2 * du);
		const double b = ddut.inner(ddur + du) + dl;
		const double s = math::sign(ddut.inner(du));
		//predictor
		return -b / a + s * sqrt(b * b - c * a) / a;
	}
	double continuation::predictor_arc_length_cylindrical(void) const
	{
		//data
		const math::vector du(m_du, m_size);
		const math::vector dur(m_dut, m_size);
		const math::vector dut(m_dur, m_size);
		//data
		const double a = dut.inner(dut);
		const double b = dut.inner(dur);
		const double s = math::sign(dut.inner(du));
		const double c = dur.inner(dur) - du.inner(du);
		//predictor
		return -b / a + s * sqrt(b * b - c * a) / a;
	}
	double continuation::corrector_arc_length_cylindrical(void) const
	{
		//data
		const math::vector du(m_du, m_size);
		const math::vector ddur(m_ddut, m_size);
		const math::vector ddut(m_ddur, m_size);
		//data
		const double a = ddut.inner(ddut);
		const double b = ddut.inner(ddur + du);
		const double c = ddur.inner(ddur + 2 * du);
		const double s = math::sign(ddut.inner(du));
		//predictor
		return -b / a + s * sqrt(b * b - c * a) / a;
	}
}