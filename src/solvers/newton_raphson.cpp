//std
#include <cmath>
#include <cstdio>
#include <cstring>

//math
#include "Math/inc/misc/misc.hpp"
#include "Math/inc/linear/matrix.hpp"
#include "Math/inc/linear/vector.hpp"
#include "Math/inc/solvers/strategies.hpp"
#include "Math/inc/solvers/newton_raphson.hpp"

namespace math
{
	//constructors
	newton_raphson::newton_raphson(void) : 
		m_silent(false), m_strategy(strategy::control_load), 
		m_stop([](){ return false; }), m_run_interface([](uint32_t){ return; }),
		m_watch_dof(0), m_step_max(100), m_attempt_max(5), m_iteration_max(10), m_nd(1), m_dl0(0.01), m_tolerance(1e-5), m_data(nullptr)
	{
		allocate();
	}

	//destructor
	newton_raphson::~newton_raphson(void)
	{
		clear();
	}

	//data
	void newton_raphson::size(uint32_t size)
	{
		clear();
		m_nd = size;
		allocate();
	}
	void newton_raphson::save(const char* path) const
	{
		FILE* file = fopen(path, "w");
		fprintf(file, "%04d\n", m_step);
		for(uint32_t i = 0; i < m_step; i++)
		{
			for(uint32_t j = 0; j < m_nd + 1; j++)
			{
				fprintf(file, "%+.6e ", m_data[j + (m_nd + 1) * i]);
			}
			fprintf(file, "\n");
		}
	}

	//solve
	void newton_raphson::clear(void)
	{
		delete[] m_data;
		m_data = nullptr;
	}
	void newton_raphson::apply(void)
	{
		m_l_new = m_l_old + m_dl;
		for(uint32_t i = 0; i < m_nd; i++)
		{
			m_x_new[i] = m_x_old[i] + m_dx[i];
		}
		m_system(m_fi.data(), m_Kt.data(), m_x_new.data());
	}
	void newton_raphson::print(void)
	{
		if(!m_silent)
		{
			printf("step: %04d load: %+.6e state: %+.6e\n", m_step, m_l_new, m_x_new[m_watch_dof]);
		}
	}
	void newton_raphson::setup(void)
	{
		//data
		m_step = 0;
		delete[] m_data;
		m_data = new double[(m_nd + 1) * m_step_max];
		m_system(m_fi.data(), m_Kt.data(), m_x_new.data());
		//state
		m_dl = m_dl0;
		m_x_old = m_x_new;
		m_l_old = m_l_new = m_fi.inner(m_fe) / m_fe.inner(m_fe);
	}
	void newton_raphson::update(void)
	{
		m_update();
		m_l_old = m_l_new;
		m_x_old = m_x_new;
		m_data[m_step * (m_nd + 1)] = m_l_new;
		for(uint32_t i = 0; i < m_nd; i++)
		{
			m_data[m_step * (m_nd + 1) + i + 1] = m_x_new[i];
		}
		m_step++;
	}
	void newton_raphson::residue(void)
	{
		//residue
		for(uint32_t i = 0; i < m_nd; i++)
		{
			m_r[i] = m_l_new * m_fe[i] - m_fi[i];
		}
		//equilibrium
		m_equilibrium = m_r.norm() < m_tolerance * m_fe.norm();
	}
	void newton_raphson::allocate(void)
	{
		m_r.resize(m_nd);
		m_dx.resize(m_nd);
		m_fi.resize(m_nd);
		m_fe.resize(m_nd);
		m_dxt.resize(m_nd);
		m_ddx.resize(m_nd);
		m_ddxr.resize(m_nd);
		m_ddxt.resize(m_nd);
		m_x_old.resize(m_nd);
		m_x_new.resize(m_nd);
		m_Kt.resize(m_nd, m_nd);
	}
	void newton_raphson::predictor(void)
	{
		m_Kt.solve(m_dxt, m_fe);
		load_predictor();
		m_dl = m_dl0 / (1 << m_attempt);
		m_dx = m_dl * m_dxt;
		apply();
	}
	void newton_raphson::corrector(void)
	{
		for(m_iteration = 0; m_iteration < m_iteration_max; m_iteration++)
		{
			//check
			residue();
			if(m_equilibrium)
			{
				break;
			}
			//corrector
			m_Kt.solve(m_ddxr, m_r);
			m_Kt.solve(m_ddxt, m_fe);
			//update
			load_corrector();
			m_dl += m_ddl;
			m_dx += m_ddxr + m_ddl * m_ddxt;
			//apply
			apply();
		}
	}
	void newton_raphson::load_predictor(void)
	{
		if(m_step != 1)
		{
			if(m_strategy == strategy::control_state)
			{
				m_dl0 = m_dx[m_watch_dof] / m_dxt[m_watch_dof];
			}
			if(m_strategy == strategy::arc_length || m_strategy == strategy::minimal_norm)
			{
				m_dl0 = sign(m_dx.inner(m_dxt)) * m_dx.norm() / m_dxt.norm();
			}
		}
	}
	void newton_raphson::load_corrector(void)
	{
		if(m_strategy == strategy::control_load)
		{
			m_ddl = 0;
		}
		if(m_strategy == strategy::control_state)
		{
			m_ddl = -m_ddxr[m_watch_dof] / m_ddxt[m_watch_dof];
		}
		if(m_strategy == strategy::arc_length)
		{
			double a = 0, b = 0, c = 0, s = 0;
			for(uint32_t i = 0; i < m_nd; i++)
			{
				s += m_ddxt[i] * m_dx[i];
				a += m_ddxt[i] * m_ddxt[i];
				b += m_ddxt[i] * (m_ddxr[i] + m_dx[i]);
				c += m_ddxr[i] * (m_ddxr[i] + 2 * m_dx[i]);
			}
			m_ddl = -b / a + sign(s) * sqrt(pow(b / a, 2) - c / a);
		}
		if(m_strategy == strategy::minimal_norm)
		{
			m_ddl = -m_ddxt.inner(m_ddxr) / m_ddxt.inner(m_ddxt);
		}
	}

	//solve
	void newton_raphson::step(void)
	{
		for(m_attempt = 0; m_attempt < m_attempt_max; m_attempt++)
		{
			predictor();
			corrector();
			if(m_equilibrium)
			{
				break;
			}
		}
		update();
	}
	void newton_raphson::solve(void)
	{
		setup();
		print();
		update();
		while(m_step < m_step_max && !m_stop())
		{
			step();
			print();
			if(!m_equilibrium)
			{
				break;
			}
			m_run_interface(m_step);
		}
	}
}