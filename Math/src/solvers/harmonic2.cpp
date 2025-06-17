//std
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <cstring>
#include <malloc.h>

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/harmonic2.hpp"

extern "C"
{
	#ifdef _WIN32
	void legendre_dr_compute(int, double[], double[]);
	#else
	void legendre_compute_dr(int, double[], double[]);
	#endif
}

#ifndef _WIN32
static void(*legendre_dr_compute)(int, double[], double[]) = legendre_compute_dr;
#endif

//x: state vector
//w: load frequency
//l: load amplitude
//r: residue vector
//v: velocity vector
//z: harmonics vector
//a: acceleration vector

//residue:
//r(t, w, l, x(t), v(t), a(t)) = l * fe(t, w, x(t), v(t)) - fi(x(t), v(t)) - M(x(t)) * a(t)

//tangent on a:
//M(t, w, l, x(t), v(t)) = M(x(t))

//tangent on v:
//C(t, w, l, x(t), v(t)) = dfi/dc(x(t), v(t)) - l * dfv/dx(t, w, x(t), v(t))

//tangent on x:
//K(t, w, l, x(t), v(t), a(t)) = dfi/dx(x(t), v(t)) - l * dfe/dx(t, w, x(t), v(t)) - dM/dx(x(t)) : a(t)

//time domain:
//z = [a0, a1, b1, ..., ak, bk]
//x(t, w, z) = a0 + cos(k * w * t) * ak + sin(k * w * t) * bk
//v(t, w, z) = -k * w * sin(k * w * t) * ak + k * w * cos(k * w * t) * bk
//a(t, w, z) = -k * k * w * w * cos(k * w * t) * ak - k * k * w * w * sin(k * w * t) * bk

//residue:
//r(t, w, l, z) = l * fe(t, w, z) - fi(t, w, z) - M(t, w, z) * a(t, w, z)

//target system:
//tk = pi / w * (1 + sk)
//u0(w, l, z) = wk * r(tk, w, l, z) = 0
//ui(w, l, z) = wk * cos(i * w * tk) * r(tk, w, l, z) = 0
//vi(w, l, z) = wk * sin(i * w * tk) * r(tk, w, l, z) = 0

namespace math
{
	namespace solvers
	{
		//constructors
		harmonic2::harmonic2(void) : 
			m_sq(nullptr), m_wq(nullptr),
			m_xd(nullptr), m_vd(nullptr), m_ad(nullptr),
			m_Kd(nullptr), m_Cd(nullptr), m_Md(nullptr)
		{
			return;
		}

		//destructor
		harmonic2::~harmonic2(void)
		{
			cleanup();
		}

		//data
		uint32_t harmonic2::state_set(void) const
		{
			return uint32_t(state::x) | uint32_t(state::p);
		}
		uint32_t harmonic2::force_set(void) const
		{
			return uint32_t(force::r) | uint32_t(force::fe);
		}
		uint32_t harmonic2::tangent_set(void) const
		{
			return uint32_t(tangent::K);
		}

		//solve
		void harmonic2::apply(void)
		{
			solver::apply();
			(m_control == control::load ? m_l : m_w) = m_p_new;
		}
		void harmonic2::check(void)
		{
			if(!m_internal_force || !m_external_force || !m_inertia || !m_damping || !m_stiffness)
			{
				printf("Error: Harmonic solver called with at least one method not set!\n");
				exit(EXIT_FAILURE);
			}
		}
		void harmonic2::setup(void)
		{
			solver::setup();
			legendre_dr_compute(m_quadrature_order, m_sq, m_wq);
		}
		void harmonic2::compute(void)
		{
			compute_harmonic_residue();
			compute_harmonic_tangent_p();
			compute_harmonic_tangent_z();
		}
		void harmonic2::predictor(void)
		{
			//data
			const matrix K(m_K, m_size, m_size);
			const vector r(m_r, m_size), fe(m_fe, m_size);
			vector dx(m_dx, m_size), dxr(m_dxr, m_size), dxt(m_dxt, m_size);
			//predictor
			if(!K.solve(dxr, r) || !K.solve(dxt, fe))
			{
				if(!m_silent) printf("Unable to decompose stiffness matrix in predictor!\n");
			}
			//continuation
			if(m_step != 1)
			{
				m_dp = m_continuation.predictor() / (1 << m_attempt);
			}
			for(uint32_t i = 0; i < m_size; i++) dx[i] = dxr[i] + m_dp * dxt[i];
		}
		void harmonic2::corrector(void)
		{
			//data
			const matrix K(m_K, m_size, m_size);
			const vector r(m_r, m_size), fe(m_fe, m_size);
			vector ddxr(m_ddxr, m_size), ddxt(m_ddxt, m_size);
			//corrector
			if(!K.solve(ddxr, r) || !K.solve(ddxt, fe))
			{
				if(!m_silent) printf("Unable to decompose stiffness matrix in corrector!\n");
			}
			//continuation
			m_ddp = m_continuation.corrector();
			//update
			m_dp += m_ddp;
			for(uint32_t i = 0; i < m_size; i++) m_dx[i] += m_ddxr[i] + m_ddp * m_ddxt[i];
		}

		//state
		void harmonic2::compute_state(double t)
		{
			memcpy(m_xd, m_x_new, m_dofs * sizeof(double));
			for(uint32_t i = 1; i <= m_harmonics; i++)
			{
				const double ci = cos(i * m_w * t);
				const double si = sin(i * m_w * t);
				const uint32_t p1 = (2 * i - 1) * m_dofs;
				const uint32_t p2 = (2 * i + 0) * m_dofs;
				for (uint32_t j = 0; j < m_dofs; j++)
				{
					m_xd[j] += ci * m_x_new[p1 + j];
					m_xd[j] += si * m_x_new[p2 + j];
				}
			}
		}
		void harmonic2::compute_residue(double t)
		{
			//state
			compute_state(t);
			compute_velocity(t);
			compute_acceleration(t);
			//state
			m_inertia(m_Md, m_xd);
			m_internal_force(m_fid, m_xd, m_vd);
			m_external_force(m_fed, m_xd, m_vd, t, m_w);
			//residue
			for(uint32_t i = 0; i < m_dofs; i++)
			{
				m_rd[i] = m_l * m_fed[i] - m_fid[i];
				for(uint32_t j = 0; j < m_dofs; j++)
				{
					m_rd[i] -= m_Md[i + m_dofs * j] * m_ad[j];
				}
			}
		}
		void harmonic2::compute_velocity(double t)
		{
			memset(m_vd, 0, m_dofs * sizeof(double));
			for(uint32_t i = 1; i <= m_harmonics; i++)
			{
				const double ci = cos(i * m_w * t);
				const double si = sin(i * m_w * t);
				const uint32_t p1 = (2 * i - 1) * m_dofs;
				const uint32_t p2 = (2 * i + 0) * m_dofs;
				for(uint32_t j = 0; j < m_dofs; j++)
				{
					m_vd[j] -= i * m_w * si * m_x_new[p1 + j];
					m_vd[j] += i * m_w * ci * m_x_new[p2 + j];
				}
			}
		}
		void harmonic2::compute_acceleration(double t)
		{
			memset(m_ad, 0, m_dofs * sizeof(double));
			for(uint32_t i = 1; i <= m_harmonics; i++)
			{
				const double ci = cos(i * m_w * t);
				const double si = sin(i * m_w * t);
				const uint32_t p1 = (2 * i - 1) * m_dofs;
				const uint32_t p2 = (2 * i + 0) * m_dofs;
				for(uint32_t j = 0; j < m_dofs; j++)
				{
					m_ad[j] -= i * i * m_w * m_w * ci * m_x_new[p1 + j];
					m_ad[j] -= i * i * m_w * m_w * si * m_x_new[p2 + j];
				}
			}
		}

		//state
		void harmonic2::compute_tangent_l(double t)
		{
			//state
			compute_state(t);
			compute_velocity(t);
			//system
			m_external_force(m_fed, m_xd, m_vd, t, m_w);
		}
		void harmonic2::compute_tangent_w(double t)
		{
			//state
			compute_state(t);
			compute_velocity(t);
			compute_acceleration(t);
			//system
			m_inertia(m_Md, m_xd);
			m_damping(m_Cd, m_xd, m_vd, t);
			memset(m_fed, 0, m_dofs * sizeof(double));
			//tangent
			for(uint32_t i = 0; i < m_dofs; i++)
			{
				for(uint32_t j = 0; j < m_dofs; j++)
				{
					m_fed[i] -= m_Cd[i + m_dofs * j] * m_vd[j] / m_w;
					m_fed[i] -= 2 * m_Md[i + m_dofs * j] * m_ad[j] / m_w;
				}
			}
		}
		void harmonic2::compute_tangent_z(double t)
		{
			//state
			compute_state(t);
			compute_velocity(t);
			compute_acceleration(t);
			//system
			m_inertia(m_Md, m_xd);
			m_damping(m_Cd, m_xd, m_vd, t);
			m_stiffness(m_Kd, m_xd, m_vd, m_ad, t, m_w, m_l);
		}

		//system
		void harmonic2::compute_harmonic_residue(void)
		{
			const uint32_t nd = m_dofs;
			const uint32_t nz = 2 * m_harmonics + 1;
			memset(m_r, 0, nd * nz * sizeof(double));
			for(uint32_t k = 0; k < m_quadrature_order; k++)
			{
				//quadrature
				const double sk = m_sq[k];
				const double wk = m_wq[k];
				const double tk = (1 + sk) * M_PI / m_w;
				//residue
				compute_residue(tk);
				for(uint32_t i = 1; i <= m_harmonics; i++)
				{
					const double ci = cos(i * m_w * tk);
					const double si = sin(i * m_w * tk);
					const uint32_t p1 = (2 * i - 1) * nd;
					const uint32_t p2 = (2 * i + 0) * nd;
					for(uint32_t j = 0; j < m_dofs; j++) m_r[p1 + j] += wk * ci * m_rd[j];
					for(uint32_t j = 0; j < m_dofs; j++) m_r[p2 + j] += wk * si * m_rd[j];
				}
				for(uint32_t i = 0; i < m_dofs; i++) m_r[i] += wk * m_rd[i];
			}
		}
		void harmonic2::compute_harmonic_tangent_p(void)
		{
			if(m_control == control::load) compute_harmonic_tangent_l();
			if(m_control == control::frequency) compute_harmonic_tangent_w();
		}
		void harmonic2::compute_harmonic_tangent_l(void)
		{
			const uint32_t nd = m_size;
			const uint32_t nz = 2 * m_harmonics + 1;
			memset(m_fe, 0, nd * nz * sizeof(double));
			for(uint32_t k = 0; k < m_quadrature_order; k++)
			{
				//quadrature
				const double sk = m_sq[k];
				const double wk = m_wq[k];
				const double tk = (1 + sk) * M_PI / m_w;
				//tangent
				compute_tangent_l(tk);
				for(uint32_t i = 1; i <= m_harmonics; i++)
				{
					const double ci = cos(i * m_w * tk);
					const double si = sin(i * m_w * tk);
					const uint32_t p1 = (2 * i - 1) * m_dofs;
					const uint32_t p2 = (2 * i + 0) * m_dofs;
					for(uint32_t j = 0; j < m_dofs; j++) m_fe[p1 + j] += wk * ci * m_fed[j];
					for(uint32_t j = 0; j < m_dofs; j++) m_fe[p2 + j] += wk * si * m_fed[j];
				}
				for(uint32_t i = 0; i < m_dofs; i++) m_fe[i] += wk * m_fed[i];
			}
		}
		void harmonic2::compute_harmonic_tangent_w(void)
		{
			const uint32_t nd = m_size;
			const uint32_t nz = 2 * m_harmonics + 1;
			memset(m_fe, 0, nd * nz * sizeof(double));
			for(uint32_t k = 0; k < m_quadrature_order; k++)
			{
				//quadrature
				const double sk = m_sq[k];
				const double wk = m_wq[k];
				const double tk = (1 + sk) * M_PI / m_w;
				//tangent
				compute_tangent_w(tk);
				for(uint32_t i = 1; i <= m_harmonics; i++)
				{
					const double ci = cos(i * m_w * tk);
					const double si = sin(i * m_w * tk);
					const uint32_t p1 = (2 * i - 1) * m_dofs;
					const uint32_t p2 = (2 * i + 0) * m_dofs;
					for(uint32_t j = 0; j < m_dofs; j++) m_fe[p1 + j] += wk * ci * m_fed[j];
					for(uint32_t j = 0; j < m_dofs; j++) m_fe[p2 + j] += wk * si * m_fed[j];
				}
				for(uint32_t j = 0; j < m_dofs; j++) m_fe[j] += wk * m_fed[j];
			}
		}
		void harmonic2::compute_harmonic_tangent_z(void)
		{
			//data
			const uint32_t nd = m_size;
			const uint32_t nz = 2 * m_harmonics + 1;
			//data
			const matrix Kd(m_Kd, nd, nd);
			const matrix Cd(m_Kd, nd, nd);
			const matrix Md(m_Kd, nd, nd);
			matrix K(m_K, nd * nz, nd * nz);
			//tangent
			K.zeros();
			for(uint32_t k = 0; k < m_quadrature_order; k++)
			{
				//quadrature
				const double sk = m_sq[k];
				const double wk = m_wq[k];
				const double tk = (1 + sk) * M_PI / m_w;
				//tangent
				compute_tangent_z(tk);
				K.span(0, 0, nd, nd) += wk * Kd;
				for(uint32_t i = 1; i <= m_harmonics; i++)
				{
					const double ci = cos(i * m_w * tk);
					const double si = sin(i * m_w * tk);
					const uint32_t p1 = (2 * i - 1) * nd;
					const uint32_t p2 = (2 * i + 0) * nd;
					K.span(p1, 0, nd, nd) += wk * ci * Kd;
					K.span(p2, 0, nd, nd) += wk * si * Kd;
					K.span(0, p1, nd, nd) -= wk * si * i * m_w * Cd;
					K.span(0, p2, nd, nd) += wk * ci * i * m_w * Cd;
					K.span(0, p1, nd, nd) += wk * ci * (Kd - i * i * m_w * m_w * Md);
					K.span(0, p2, nd, nd) += wk * si * (Kd - i * i * m_w * m_w * Md);
					for(uint32_t j = 1; j <= m_harmonics; j++)
					{
						const double cj = cos(j * m_w * tk);
						const double sj = sin(j * m_w * tk);
						const uint32_t q1 = (2 * j - 1) * nd;
						const uint32_t q2 = (2 * j + 0) * nd;
						K.span(p1, q1, nd, nd) -= wk * ci * sj * j * m_w * Cd;
						K.span(p1, q2, nd, nd) += wk * ci * cj * j * m_w * Cd;
						K.span(p2, q1, nd, nd) -= wk * si * sj * j * m_w * Cd;
						K.span(p2, q2, nd, nd) += wk * si * cj * j * m_w * Cd;
						K.span(p1, q1, nd, nd) += wk * ci * cj * (Kd - j * j * m_w * m_w * Md);
						K.span(p1, q2, nd, nd) += wk * ci * sj * (Kd - j * j * m_w * m_w * Md);
						K.span(p2, q1, nd, nd) += wk * si * cj * (Kd - j * j * m_w * m_w * Md);
						K.span(p2, q2, nd, nd) += wk * si * sj * (Kd - j * j * m_w * m_w * Md);
					}
				}
			}
		}

		//solver
		void harmonic2::cleanup(void)
		{
			//harmonic
			double* data[] = {
				m_sq, m_wq, m_xd, m_vd, m_ad, m_rd, m_fid, m_fed, m_Kd, m_Cd, m_Md
			};
			for(double* ptr : data)
			{
				delete[] ptr;
			}
			//solver
			solver::cleanup();
		}
		void harmonic2::allocate(void)
		{
			//harmonic
			m_xd = new double[m_dofs];
			m_vd = new double[m_dofs];
			m_ad = new double[m_dofs];
			m_rd = new double[m_dofs];
			m_fid = new double[m_dofs];
			m_fed = new double[m_dofs];
			m_Kd = new double[m_dofs * m_dofs];
			m_Cd = new double[m_dofs * m_dofs];
			m_Md = new double[m_dofs * m_dofs];
			m_sq = new double[m_quadrature_order];
			m_wq = new double[m_quadrature_order];
			m_size = (1 + 2 * m_harmonics) * m_dofs;
			//solver
			solver::allocate();
		}
	}
}