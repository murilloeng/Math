//std
#include <cstdio>
#include <stdexcept>

//Math
#include "Math/inc/Solvers/NewtonRaphson.hpp"

//x: state Vector
//r: residual force Vector
//p: continuation parameter

//target system: r(x, p) = 0

//tangent on x: K(x, p) = -dr/dx(x, p)
//tangent on p: g(x, p) = +dr/dp(x, p)

namespace math
{
	namespace solvers
	{
		//constructors
		NewtonRaphson::NewtonRaphson(void)
		{
			return;
		}

		//destructor
		NewtonRaphson::~NewtonRaphson(void)
		{
			return;
		}

		//solve
		void NewtonRaphson::compute_step(void)
		{
			for(m_attempt = 0; m_attempt < m_attempt_max; m_attempt++)
			{
				predictor();
				for(m_iteration = 0; m_iteration < m_iteration_max; m_iteration++)
				{
					apply();
					compute();
					if(equilibrium()) break; else corrector();
				}
				if(m_status) break;
				restore();
			}
			update();
			record();
		}

		//data
		uint32_t NewtonRaphson::state_set(void) const
		{
			return 1 << uint32_t(State::x) | 1 << uint32_t(State::p);
		}
		uint32_t NewtonRaphson::force_set(void) const
		{
			return 1 << uint32_t(Force::r) | 1 << uint32_t(Force::fe);
		}
		uint32_t NewtonRaphson::tangent_set(void) const
		{
			return 1 << uint32_t(Tangent::K);
		}

		//data
		NewtonRaphson::System NewtonRaphson::system(void) const
		{
			return m_system;
		}
		NewtonRaphson::System NewtonRaphson::system(NewtonRaphson::System system)
		{
			return m_system = system;
		}

		NewtonRaphson::Function NewtonRaphson::residue(void) const
		{
			return m_residue;
		}
		NewtonRaphson::Function NewtonRaphson::residue(NewtonRaphson::Function residue)
		{
			return m_residue = residue;
		}

		NewtonRaphson::Function NewtonRaphson::tangent_p(void) const
		{
			return m_tangent_p;
		}
		NewtonRaphson::Function NewtonRaphson::tangent_p(NewtonRaphson::Function tangent_p)
		{
			return m_tangent_p = tangent_p;
		}

		NewtonRaphson::Function NewtonRaphson::tangent_x(void) const
		{
			return m_tangent_x;
		}
		NewtonRaphson::Function NewtonRaphson::tangent_x(NewtonRaphson::Function tangent_x)
		{
			return m_tangent_x = tangent_x;
		}

		//solve
		void NewtonRaphson::check(void)
		{
			if(!m_system && !(m_residue && m_tangent_p && m_tangent_x))
			{
				throw std::runtime_error("Newton-Raphson solver called with at least one method not set!");
			}
		}
		void NewtonRaphson::print(void)
		{
			Incremental::print();
			Implicit::print();
			Solver::print();
		}
		void NewtonRaphson::setup(void)
		{
			Solver::setup();
			Implicit::setup();
			Incremental::setup();
		}
		void NewtonRaphson::compute(void)
		{
			if(m_system)
			{
				m_system(m_r, m_fe, m_K, m_p_new, m_x_new);
			}
			else
			{
				m_residue(m_r, m_p_new, m_x_new);
				m_tangent_x(m_K, m_p_new, m_x_new);
				m_tangent_p(m_fe, m_p_new, m_x_new);
			}
		}
		void NewtonRaphson::predictor(void)
		{
			//predictor
			if(!Solver::solve(m_K, m_r, m_dxr) || !Solver::solve(m_K, m_fe, m_dxt))
			{
				if(!m_silent) printf("Unable to decompose stiffness Matrix in predictor!\n");
			}
			//continuation
			if(m_step != 1)
			{
				m_dp = m_continuation.predictor() / (1 << m_attempt);
			}
			for(uint32_t i = 0; i < m_size; i++) m_dx[i] = m_dxr[i] + m_dp * m_dxt[i];
		}
		void NewtonRaphson::corrector(void)
		{
			//corrector
			if(!Solver::solve(m_K, m_r, m_ddxr) || !Solver::solve(m_K, m_fe, m_ddxt))
			{
				if(!m_silent) printf("Unable to decompose stiffness Matrix in corrector!\n");
			}
			m_ddp = m_continuation.corrector();
			//update
			m_dp += m_ddp;
			for(uint32_t i = 0; i < m_size; i++) m_dx[i] += m_ddxr[i] + m_ddp * m_ddxt[i];
		}
	}
}