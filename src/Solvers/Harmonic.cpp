//std
#include <cmath>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Miscellaneous/util.hpp"
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Solvers/Harmonic.hpp"

extern "C"
{
	void legendre_compute_dr(int, double[], double[]);
}

//x: state Vector
//w: load frequency
//l: load amplitude
//r: residue Vector
//v: velocity Vector
//z: harmonics Vector
//a: acceleration Vector

//residue:
//r(t, w, l, x(t), v(t), a(t)) = l * fe(t, w, x(t), v(t)) - fi(x(t), v(t)) - M(x(t)) * a(t)

//tangent on a:
//M(t, w, l, x(t), v(t)) = M(x(t))

//tangent on v:
//C(t, w, l, x(t), v(t)) = dfi/dv(x(t), v(t)) - l * dfe/dv(t, w, x(t), v(t))

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
		Harmonic::Harmonic(void) : 
			m_sq{nullptr}, m_wq{nullptr},
			m_xd{nullptr}, m_vd{nullptr}, m_ad{nullptr},
			m_Kd{nullptr}, m_Cd{nullptr}, m_Md{nullptr},
			m_rd{nullptr}, m_fid{nullptr}, m_fed{nullptr},
			m_dofs{0}, m_harmonics{0}, m_quadrature_order{20},
			m_stability_steps{100}, m_stability{false}, m_stability_data{nullptr}
		{
			return;
		}

		//destructor
		Harmonic::~Harmonic(void)
		{
			const double* data[] = {
				m_sq, m_wq, m_xd, m_vd, m_ad, m_Kd, m_Cd, m_Md, m_rd, m_fid, m_fed
			};
			for(const double* ptr : data)
			{
				delete[] ptr;
			}
			delete[] m_stability_data;
		}

		//tests
		void Harmonic::test_inertia(void) const
		{
			//data
			double t, w, l;
			const uint32_t nt = 1000;
			const uint32_t nd = m_dofs;
			math::Vector xd(nd), vd(nd), ad(nd);
			math::Matrix Ma(nd, nd), Mn(nd, nd), Me(nd, nd);
			//function
			std::function<void(double*, const double*)> function = [this, &t, &w, &l, nd, &xd, &vd](double* r, const double* a){
				math::Matrix M(nd, nd);
				math::Vector fi(nd), fe(nd);
				m_inertia(M.data(), xd.data());
				m_external_force(fe.data(), xd.data(), t, w);
				m_internal_force(fi.data(), xd.data(), vd.data());
				math::Vector(r, nd) = l * fe - fi - M * math::Vector(a, nd);
			};
			//gradient
			std::function<void(double*, const double*)> gradient = [this, nd, &xd](double* M, const double* a){
				m_inertia(M, xd.data());
				math::Matrix(M, nd, nd) *= -1;
			};
			//test
			for(uint32_t i = 0; i < nt; i++)
			{
				//state
				xd.randu();
				vd.randu();
				ad.randu();
				t = math::randu(0, 1);
				w = math::randu(0, 1);
				l = math::randu(0, 1);
				//test
				gradient(Ma.data(), ad.data());
				math::ndiff(function, Mn.data(), ad.data(), nd, nd, 1e-5);
				//check
				Me = Ma - Mn;
				if(Me.norm() > 1e-5 * Ma.norm())
				{
					Ma.print("Ma");
					Mn.print("Mn");
					Me.print("Me");
					break;
				}
				//print
				printf("test: %d error norm: %+.2e\n", i, Me.norm());
			}
		}
		void Harmonic::test_damping(void) const
		{
			//data
			double t, w, l;
			const uint32_t nt = 1000;
			const uint32_t nd = m_dofs;
			math::Vector xd(nd), vd(nd), ad(nd);
			math::Matrix Ca(nd, nd), Cn(nd, nd), Ce(nd, nd);
			//function
			std::function<void(double*, const double*)> function = [this, &t, &w, &l, nd, &xd, &ad](double* r, const double* v){
				math::Matrix M(nd, nd);
				math::Vector fi(nd), fe(nd);
				m_inertia(M.data(), xd.data());
				m_internal_force(fi.data(), xd.data(), v);
				m_external_force(fe.data(), xd.data(), t, w);
				math::Vector(r, nd) = l * fe - fi - M * ad;
			};
			//gradient
			std::function<void(double*, const double*)> gradient = [this, nd, &xd](double* C, const double* v){
				m_damping(C, xd.data(), v);
				math::Matrix(C, nd, nd) *= -1;
			};
			//test
			for(uint32_t i = 0; i < nt; i++)
			{
				//state
				xd.randu();
				vd.randu();
				ad.randu();
				t = math::randu(0, 1);
				w = math::randu(0, 1);
				l = math::randu(0, 1);
				//test
				gradient(Ca.data(), vd.data());
				math::ndiff(function, Cn.data(), vd.data(), nd, nd, 1e-5);
				//check
				Ce = Ca - Cn;
				if(Ce.norm() > 1e-5 * Ca.norm())
				{
					Ca.print("Ca");
					Cn.print("Cn");
					Ce.print("Ce");
					break;
				}
				//print
				printf("test: %d error norm: %+.2e\n", i, Ce.norm());
			}
		}
		void Harmonic::test_stiffness(void) const
		{
			//data
			double t, w, l;
			const uint32_t nt = 1000;
			const uint32_t nd = m_dofs;
			math::Vector xd(nd), vd(nd), ad(nd);
			math::Matrix Ka(nd, nd), Kn(nd, nd), Ke(nd, nd);
			//function
			std::function<void(double*, const double*)> function = [this, &t, &w, &l, nd, &vd, &ad](double* r, const double* x){
				math::Matrix M(nd, nd);
				math::Vector fi(nd), fe(nd);
				m_inertia(M.data(), x);
				m_external_force(fe.data(), x, t, w);
				m_internal_force(fi.data(), x, vd.data());
				math::Vector(r, nd) = l * fe - fi - M * ad;
			};
			//gradient
			std::function<void(double*, const double*)> gradient = [this, nd, &vd, &ad, &t, &w, &l](double* K, const double* x){
				m_stiffness(K, x, vd.data(), ad.data(), t, w, l);
				math::Matrix(K, nd, nd) *= -1;
			};
			//test
			for(uint32_t i = 0; i < nt; i++)
			{
				//state
				xd.randu();
				vd.randu();
				ad.randu();
				t = math::randu(0, 1);
				w = math::randu(0, 1);
				l = math::randu(0, 1);
				//test
				gradient(Ka.data(), vd.data());
				math::ndiff(function, Kn.data(), vd.data(), nd, nd, 1e-5);
				//check
				Ke = Ka - Kn;
				if(Ke.norm() > 1e-5 * Ka.norm())
				{
					Ka.print("Ka");
					Kn.print("Kn");
					Ke.print("Ke");
					break;
				}
				//print
				printf("test: %d error norm: %+.2e\n", i, Ke.norm());
			}
		}

		//solve
		void Harmonic::apply(void)
		{
			Solver::apply();
			(m_control == Control::Load ? m_l : m_w) = m_p_new;
		}
		void Harmonic::check(void)
		{
			if(!m_internal_force || !m_external_force || !m_inertia || !m_damping || !m_stiffness)
			{
				throw std::runtime_error("Error: Harmonic solver called with at least one method not set!");
			}
		}
		void Harmonic::setup(void)
		{
			//data
			m_p_new = m_control == Control::Load ? m_l : m_w;
			legendre_compute_dr(m_quadrature_order, m_sq, m_wq);
			//initial
			m_dp = 0;
			m_convergence.m_solver = this;
			for(m_iteration = 0; m_iteration < m_iteration_max; m_iteration++)
			{
				//data
				compute();
				if(equilibrium()) break;
				//corrector
				if(!solve(m_K, m_r, m_dx))
				{
					if(!m_silent) printf("Unable to decompose stiffness Matrix in setup!\n");
				}
				for(uint32_t i = 0; i < m_size; i++) m_x_new[i] += m_dx[i];
			}
			//setup
			Solver::setup();
		}

		//state
		void Harmonic::compute_state(const double* z, double t)
		{
			memcpy(m_xd, z, m_dofs * sizeof(double));
			for(uint32_t i = 1; i <= m_harmonics; i++)
			{
				const double ci = cos(i * m_w * t);
				const double si = sin(i * m_w * t);
				const uint32_t p1 = (2 * i - 1) * m_dofs;
				const uint32_t p2 = (2 * i + 0) * m_dofs;
				for(uint32_t j = 0; j < m_dofs; j++)
				{
					m_xd[j] += ci * z[p1 + j];
					m_xd[j] += si * z[p2 + j];
				}
			}
		}
		void Harmonic::compute_residue(const double* z, double t)
		{
			//state
			compute_state(z, t);
			compute_velocity(z, t);
			compute_acceleration(z, t);
			//state
			m_inertia(m_Md, m_xd);
			m_internal_force(m_fid, m_xd, m_vd);
			m_external_force(m_fed, m_xd, t, m_w);
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
		void Harmonic::compute_velocity(const double* z, double t)
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
					m_vd[j] -= i * m_w * si * z[p1 + j];
					m_vd[j] += i * m_w * ci * z[p2 + j];
				}
			}
		}
		void Harmonic::compute_acceleration(const double* z, double t)
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
					m_ad[j] -= i * i * m_w * m_w * ci * z[p1 + j];
					m_ad[j] -= i * i * m_w * m_w * si * z[p2 + j];
				}
			}
		}

		//state
		void Harmonic::compute_tangent_l(const double* z, double t)
		{
			//state
			compute_state(z, t);
			compute_velocity(z, t);
			//system
			m_external_force(m_fed, m_xd, t, m_w);
		}
		void Harmonic::compute_tangent_w(const double* z, double t)
		{
			//state
			compute_state(z, t);
			compute_velocity(z, t);
			compute_acceleration(z, t);
			//system
			m_inertia(m_Md, m_xd);
			m_damping(m_Cd, m_xd, m_vd);
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
		void Harmonic::compute_tangent_z(const double* z, double t)
		{
			//state
			compute_state(z, t);
			compute_velocity(z, t);
			compute_acceleration(z, t);
			//system
			m_inertia(m_Md, m_xd);
			m_damping(m_Cd, m_xd, m_vd);
			m_stiffness(m_Kd, m_xd, m_vd, m_ad, t, m_w, m_l);
		}

		//system
		void Harmonic::compute_harmonic_residue(double* r, const double* z)
		{
			memset(r, 0, m_size * sizeof(double));
			for(uint32_t k = 0; k < m_quadrature_order; k++)
			{
				//quadrature
				const double sk = m_sq[k];
				const double wk = m_wq[k];
				const double tk = (1 + sk) * M_PI / m_w;
				//residue
				compute_residue(z, tk);
				for(uint32_t i = 1; i <= m_harmonics; i++)
				{
					const double ci = cos(i * m_w * tk);
					const double si = sin(i * m_w * tk);
					const uint32_t p1 = (2 * i - 1) * m_dofs;
					const uint32_t p2 = (2 * i + 0) * m_dofs;
					for(uint32_t j = 0; j < m_dofs; j++) r[p1 + j] += wk * ci * m_rd[j];
					for(uint32_t j = 0; j < m_dofs; j++) r[p2 + j] += wk * si * m_rd[j];
				}
				for(uint32_t i = 0; i < m_dofs; i++) r[i] += wk * m_rd[i];
			}
		}
		void Harmonic::compute_harmonic_tangent_p(double* g, const double* d)
		{
			if(m_control == Control::Load) compute_harmonic_tangent_l(g, d);
			if(m_control == Control::Frequency) compute_harmonic_tangent_w(g, d);
		}
		void Harmonic::compute_harmonic_tangent_l(double* g, const double* z)
		{
			memset(g, 0, m_size * sizeof(double));
			for(uint32_t k = 0; k < m_quadrature_order; k++)
			{
				//quadrature
				const double sk = m_sq[k];
				const double wk = m_wq[k];
				const double tk = (1 + sk) * M_PI / m_w;
				//tangent
				compute_tangent_l(z, tk);
				for(uint32_t i = 1; i <= m_harmonics; i++)
				{
					const double ci = cos(i * m_w * tk);
					const double si = sin(i * m_w * tk);
					const uint32_t p1 = (2 * i - 1) * m_dofs;
					const uint32_t p2 = (2 * i + 0) * m_dofs;
					for(uint32_t j = 0; j < m_dofs; j++) g[p1 + j] += wk * ci * m_fed[j];
					for(uint32_t j = 0; j < m_dofs; j++) g[p2 + j] += wk * si * m_fed[j];
				}
				for(uint32_t i = 0; i < m_dofs; i++) g[i] += wk * m_fed[i];
			}
		}
		void Harmonic::compute_harmonic_tangent_w(double* g, const double* z)
		{
			memset(g, 0, m_size * sizeof(double));
			for(uint32_t k = 0; k < m_quadrature_order; k++)
			{
				//quadrature
				const double sk = m_sq[k];
				const double wk = m_wq[k];
				const double tk = (1 + sk) * M_PI / m_w;
				//tangent
				compute_tangent_w(z, tk);
				for(uint32_t i = 1; i <= m_harmonics; i++)
				{
					const double ci = cos(i * m_w * tk);
					const double si = sin(i * m_w * tk);
					const uint32_t p1 = (2 * i - 1) * m_dofs;
					const uint32_t p2 = (2 * i + 0) * m_dofs;
					for(uint32_t j = 0; j < m_dofs; j++) g[p1 + j] += wk * ci * m_fed[j];
					for(uint32_t j = 0; j < m_dofs; j++) g[p2 + j] += wk * si * m_fed[j];
				}
				for(uint32_t i = 0; i < m_dofs; i++) g[i] += wk * m_fed[i];
			}
		}
		void Harmonic::compute_harmonic_tangent_z(double* K, const double* z)
		{
			//data
			Matrix Ks(K, m_size, m_size);
			const Matrix Kd(m_Kd, m_dofs, m_dofs);
			const Matrix Cd(m_Cd, m_dofs, m_dofs);
			const Matrix Md(m_Md, m_dofs, m_dofs);
			//tangent
			memset(K, 0, m_size * m_size * sizeof(double));
			for(uint32_t k = 0; k < m_quadrature_order; k++)
			{
				//quadrature
				const double sk = m_sq[k];
				const double wk = m_wq[k];
				const double tk = (1 + sk) * M_PI / m_w;
				//tangent
				compute_tangent_z(z, tk);
				Ks.span(0, 0, m_dofs, m_dofs) += wk * Kd;
				for(uint32_t i = 1; i <= m_harmonics; i++)
				{
					const double ci = cos(i * m_w * tk);
					const double si = sin(i * m_w * tk);
					const uint32_t p1 = (2 * i - 1) * m_dofs;
					const uint32_t p2 = (2 * i + 0) * m_dofs;
					Ks.span(p1, 0, m_dofs, m_dofs) += wk * ci * Kd;
					Ks.span(p2, 0, m_dofs, m_dofs) += wk * si * Kd;
					Ks.span(0, p1, m_dofs, m_dofs) -= wk * si * i * m_w * Cd;
					Ks.span(0, p2, m_dofs, m_dofs) += wk * ci * i * m_w * Cd;
					Ks.span(0, p1, m_dofs, m_dofs) += wk * ci * (Kd - i * i * m_w * m_w * Md);
					Ks.span(0, p2, m_dofs, m_dofs) += wk * si * (Kd - i * i * m_w * m_w * Md);
					for(uint32_t j = 1; j <= m_harmonics; j++)
					{
						const double cj = cos(j * m_w * tk);
						const double sj = sin(j * m_w * tk);
						const uint32_t q1 = (2 * j - 1) * m_dofs;
						const uint32_t q2 = (2 * j + 0) * m_dofs;
						Ks.span(p1, q1, m_dofs, m_dofs) -= wk * ci * sj * j * m_w * Cd;
						Ks.span(p1, q2, m_dofs, m_dofs) += wk * ci * cj * j * m_w * Cd;
						Ks.span(p2, q1, m_dofs, m_dofs) -= wk * si * sj * j * m_w * Cd;
						Ks.span(p2, q2, m_dofs, m_dofs) += wk * si * cj * j * m_w * Cd;
						Ks.span(p1, q1, m_dofs, m_dofs) += wk * ci * cj * (Kd - j * j * m_w * m_w * Md);
						Ks.span(p1, q2, m_dofs, m_dofs) += wk * ci * sj * (Kd - j * j * m_w * m_w * Md);
						Ks.span(p2, q1, m_dofs, m_dofs) += wk * si * cj * (Kd - j * j * m_w * m_w * Md);
						Ks.span(p2, q2, m_dofs, m_dofs) += wk * si * sj * (Kd - j * j * m_w * m_w * Md);
					}
				}
			}
		}

		//solver
		void Harmonic::solve(void)
		{
			//setup
			m_system_2 = [this](double* r, double* g, double* K, double p, const double* z){
				compute_harmonic_residue(r, z);
				compute_harmonic_tangent_p(g, z);
				compute_harmonic_tangent_z(K, z);
			};
			//solve
			NewtonRaphson::solve();
		}
		void Harmonic::cleanup(void)
		{
			//harmonic
			double** data[] = {
				&m_sq, &m_wq, &m_xd, &m_vd, &m_ad, &m_Kd, &m_Cd, &m_Md, &m_rd, &m_fid, &m_fed
			};
			for(double** ptr : data)
			{
				delete[] *ptr;
				*ptr = nullptr;
			}
			//stability
			delete[] m_stability_data;
			m_stability_data = nullptr;
			//solver
			NewtonRaphson::cleanup();
		}
		void Harmonic::allocate(void)
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
			if(m_stability) m_stability_data = new bool[m_step_max + 1];
			//solver
			NewtonRaphson::allocate();
			memset(m_x_new, 0, m_size * sizeof(double));
			memset(m_x_old, 0, m_size * sizeof(double));
		}
	}
}