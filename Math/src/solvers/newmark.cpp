//std
#include <cstdio>
#include <cstring>

//math
#include "Math/Math/inc/solvers/newmark.hpp"

namespace math
{
	//constructors
	newmark::newmark(uint32_t m_nd, bool mem) : m_mem(mem), m_nd(m_nd), m_g(0.50), m_b(0.25)
	{
		//state
		m_x = m_mem ? new double[m_nd] : nullptr;
		m_v = m_mem ? new double[m_nd] : nullptr;
		m_a = m_mem ? new double[m_nd] : nullptr;
		//force
		m_r = m_mem ? new double[m_nd] : nullptr;
		m_fi = m_mem ? new double[m_nd] : nullptr;
		m_fe = m_mem ? new double[m_nd] : nullptr;
		//increment
		m_dx = m_mem ? new double[m_nd] : nullptr;
		m_dv = m_mem ? new double[m_nd] : nullptr;
		//tangent
		m_K = m_mem ? new double[m_nd * m_nd] : nullptr;
		m_C = m_mem ? new double[m_nd * m_nd] : nullptr;
		m_M = m_mem ? new double[m_nd * m_nd] : nullptr;
	}

	//destructor
	newmark::~newmark(void)
	{
		if(m_mem)
		{
			delete[] m_x;
			delete[] m_v;
			delete[] m_a;
			delete[] m_r;
			delete[] m_K;
			delete[] m_C;
			delete[] m_M;
			delete[] m_fi;
			delete[] m_fe;
			delete[] m_dx;
			delete[] m_dv;
		}
	}

	//data
	void newmark::update(void)
	{
		inertia();
		internal();
		external();
		for(uint32_t i = 0; i < m_nd; i++)
		{
			m_r[i] = m_fe[i] - m_fi[i];
		}
		// mat::solve(m_a, m_M, m_r, m_nd);
	}
	void newmark::residue(void)
	{
		inertia();
		internal();
		external();
		for(uint32_t i = 0; i < m_nd; i++)
		{
			m_r[i] = m_fe[i] - m_fi[i];
			for(uint32_t j = 0; j < m_nd; j++)
			{
				m_r[i] -= m_M[i + m_nd * j] * m_a[j];
			}
		}
	}
	void newmark::inertia(void)
	{
		m_inertia(m_M, m_x);
	}
	void newmark::damping(void)
	{
		m_damping(m_C, m_x, m_v);
	}
	void newmark::stifness(void)
	{
		m_stiffness(m_K, m_x, m_v);
	}
	void newmark::internal(void)
	{
		m_internal(m_fi, m_x, m_v);
	}
	void newmark::external(void)
	{
		m_external(m_fe, m_t);
	}

	//solve
	void newmark::setup(void)
	{
		m_t = 0;
		m_s = 0;
		m_dt = m_T / m_ns;
	}
	void newmark::predictor(void)
	{
		m_t += m_dt;
		for(uint32_t i = 0; i < m_nd; i++)
		{
			m_v[i] += m_dt * m_a[i];
			m_x[i] += m_dt * m_v[i] - m_dt * m_dt / 2 * m_a[i];
		}
	}
	void newmark::corrector(void)
	{
		while(true)
		{
			residue();
			// const double f = mat::norm(m_fi, m_nd);
			// if(mat::norm(m_r, m_nd) < 1e-5 * (f == 0 ? 1 : f))
			// {
			// 	break;
			// }
			damping();
			stifness();
			for(uint32_t i = 0; i < m_nd * m_nd; i++)
			{
				m_K[i] += (m_g * m_dt * m_C[i] + m_M[i]) / (m_b * m_dt * m_dt);
			}
			// mat::solve(m_dx, m_K, m_r, m_nd);
			for(uint32_t i = 0; i < m_nd; i++)
			{
				m_x[i] += m_dx[i];
				m_v[i] += m_dx[i] * m_g / (m_b * m_dt);
				m_a[i] += m_dx[i] / (m_b * m_dt * m_dt);
			}
		}
	}
	void newmark::serialize(void)
	{
		printf("%04d ", m_s);
		printf("%+.6e ", m_t);
		for(uint32_t i = 0; i < m_nd; i++)
		{
			printf("%+.6e %+.6e %+.6e ", m_x[i], m_v[i], m_a[i]);
		}
		printf("\n");
		m_s++;
	}

	//solve
	void newmark::step(void)
	{
		predictor();
		corrector();
	}
	void newmark::solve(void)
	{
		setup();
		update();
		serialize();
		while(m_s < m_ns)
		{
			step();
			serialize();
		}
	}
}