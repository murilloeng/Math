//std
#include <cstdio>
#include <cstring>

//Math
#include "Math/inc/solvers/RungeKutta.hpp"

namespace math
{
	namespace solvers
	{
		//constructors
		RungeKutta::RungeKutta(uint32_t nd, bool mem) : m_mem(mem), m_type(true), m_nd(nd), m_ns(100), m_T(100)
		{
			m_x = m_mem ? new double[m_nd] : nullptr;
			m_v = m_mem ? new double[m_nd] : nullptr;
			m_a = m_mem ? new double[m_nd] : nullptr;
			m_dx = m_mem ? new double[m_nd] : nullptr;
			m_dv = m_mem ? new double[m_nd] : nullptr;
			m_xn = m_mem ? new double[m_nd] : nullptr;
			m_vn = m_mem ? new double[m_nd] : nullptr;
		}

		//destructor
		RungeKutta::~RungeKutta(void)
		{
			if(m_mem)
			{
				delete[] m_x;
				delete[] m_v;
				delete[] m_a;
				delete[] m_xn;
				delete[] m_vn;
				delete[] m_dx;
				delete[] m_dv;
			}
		}

		//data
		void RungeKutta::update(void)
		{
			!m_type ? m_system_1(m_v, m_x, m_t) : m_system_2(m_a, m_x, m_v, m_t);
		}
		void RungeKutta::state(double f)
		{
			if(!m_type)
			{
				for(uint32_t i = 0; i < m_nd; i++)
				{
					m_x[i] = m_xn[i] + f * m_v[i];
				}
			}
			else
			{
				for(uint32_t i = 0; i < m_nd; i++)
				{
					m_x[i] = m_xn[i] + f * m_v[i];
					m_v[i] = m_vn[i] + f * m_a[i];
				}
			}
		}
		void RungeKutta::increment(double f)
		{
			if(!m_type)
			{
				for(uint32_t i = 0; i < m_nd; i++)
				{
					m_dx[i] += f * m_v[i];
				}
			}
			else
			{
				for(uint32_t i = 0; i < m_nd; i++)
				{
					m_dx[i] += f * m_v[i];
					m_dv[i] += f * m_a[i];
				}
			}
		}

		//solve
		void RungeKutta::setup(void)
		{
			m_t = 0;
			m_s = 0;
			m_dt = m_T / m_ns;
			memcpy(m_xn, m_x, m_nd * sizeof(double));
			memcpy(m_vn, m_v, m_nd * sizeof(double));
		}
		void RungeKutta::tangent_1(void)
		{
			memset(m_dx, 0, m_nd * sizeof(double));
			memset(m_dv, 0, m_nd * sizeof(double));
			increment(m_dt / 6);
		}
		void RungeKutta::tangent_2(void)
		{
			//time
			m_t += m_dt / 2;
			//update state
			state(m_dt / 2);
			//update tangent
			update();
			increment(m_dt / 3);
		}
		void RungeKutta::tangent_3(void)
		{
			//update state
			state(m_dt / 2);
			//update tangent
			update();
			increment(m_dt / 3);
		}
		void RungeKutta::tangent_4(void)
		{
			//time
			m_t += m_dt / 2;
			//update state
			state(m_dt);
			//update tangent
			update();
			increment(m_dt / 6);
		}
		void RungeKutta::corrector(void)
		{
			//update state
			for(uint32_t i = 0; i < m_nd; i++)
			{
				m_x[i] = m_xn[i] += m_dx[i];
				m_v[i] = m_vn[i] += m_dv[i];
			}
			//update system
			update();
		}

		//solve
		void RungeKutta::init(void)
		{
			setup();
			update();
			serialize();
		}
		void RungeKutta::step(void)
		{
			tangent_1();
			tangent_2();
			tangent_3();
			tangent_4();
			corrector();
		}
		void RungeKutta::solve(void)
		{
			init();
			while(m_s < m_ns)
			{
				step();
				serialize();
			}
		}
		void RungeKutta::serialize(void)
		{
			m_s++;
			printf("%04d %+.6e ", m_s, m_t);
			for(uint32_t i = 0; i < m_nd; i++)
			{
				printf("%+.6e %+.6e %+.6e ", m_x[i], m_v[i], m_a[i]);
			}
			printf("\n");
		}
	}
}