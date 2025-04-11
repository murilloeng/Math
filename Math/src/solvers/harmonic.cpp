//std
#include <cmath>
#include <cstdio>
#include <cstring>

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/harmonic.hpp"

extern "C"
{
	void legendre_dr_compute(int, double[], double[]);
}

namespace math
{
	//constructors
	harmonic::harmonic(void) : 
		m_args(nullptr), m_tolerance(0),
		m_size(0), m_step_max(0), m_harmonics(0), m_attempt_max(0), m_iteration_max(0), 
		m_quadrature_order(0), m_parameter(harmonic_parameter::load),
		m_internal_force(nullptr), m_external_force(nullptr), m_external_force_gradient(nullptr),
		m_inertia(nullptr), m_damping(nullptr), m_stiffness(nullptr),
		m_sq(nullptr), m_wq(nullptr), 
		m_d(nullptr), m_v(nullptr), m_a(nullptr), m_j(nullptr),
		m_ddw(nullptr), m_dvw(nullptr), m_daw(nullptr),
		m_l_old(0), m_l_new(0), m_l_data(nullptr),
		m_w_old(0), m_w_new(0), m_w_data(nullptr),
		m_z_old(nullptr), m_z_new(nullptr), m_z_data(nullptr),
		m_dp(0), m_dpg(0), m_ddp(0),
		m_dz(nullptr), m_dz0r(nullptr), m_dz0t(nullptr), m_ddzr(nullptr), m_ddzt(nullptr),
		m_r(nullptr), m_fi(nullptr), m_fe(nullptr), m_fr(nullptr),
		m_Kt(nullptr), m_Ct(nullptr), m_Mt(nullptr), m_At(nullptr), m_bt(nullptr), m_dfew(nullptr), m_dfrw(nullptr)
	{
		return;
	}

	//destructor
	harmonic::~harmonic(void)
	{
		cleanup();
	}

	//solve
	bool harmonic::solve(void)
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//solve
		setup();
		for(m_step = 0; m_step < m_step_max; m_step++)
		{
			compute_predictor();
			for(m_iteration = 0; m_iteration < m_iteration_max; m_iteration++)
			{
				//equilibrium
				compute_system_residue();
				m_system_equilibrium = vector(m_r, nd * nz).norm() < m_tolerance;
				//check
				if(m_system_equilibrium) break;
				//corrector
				compute_corrector();
			}
			if(!m_system_equilibrium)
			{
				printf("Harmonic solver failed to converge at step %d\n", m_step);
				return false;
			}
			else
			{
				record();
				update();
				printf("step: %04d : iteration: %02d\n", m_step, m_iteration);
			}
		}
		//return
		return true;
	}

	//data
	const double* harmonic::amplitudes(void) const
	{
		return m_z_new;
	}

	//solver
	void harmonic::setup(void)
	{
		//data
		cleanup();
		allocate();
		//loop
		m_step = 0;
		m_attempt = 0;
		m_iteration = 0;
		//quadrature
		legendre_dr_compute(m_quadrature_order, m_sq, m_wq);
	}
	void harmonic::cleanup(void)
	{
		double* data[] = {
			m_sq, m_wq,
			m_d, m_v, m_a, m_j, m_ddw, m_dvw, m_daw,
			m_l_data, m_w_data, m_z_old, m_z_new, m_z_data,
			m_dz, m_dz0r, m_dz0t, m_ddzr, m_ddzt,
			m_r, m_fi, m_fe, m_fr, m_Kt, m_Ct, m_Mt, m_At, m_bt, m_dfew, m_dfrw
		};
		for(double* ptr : data) delete[] ptr;
	}
	void harmonic::allocate(void)
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t ns = m_step_max;
		const uint32_t nq = m_quadrature_order;
		const uint32_t nz = 2 * m_harmonics + 1;
		//allocate nd
		double** data_nd[] = {
			&m_d, &m_v, &m_a, &m_j, &m_fi, &m_fe, &m_fr, &m_ddw, &m_dvw, &m_daw, &m_dfew, &m_dfrw
		};
		for(double** ptr : data_nd) *ptr = new double[nd];
		//allocate nq
		m_sq = new double[nq];
		m_wq = new double[nq];
		//allocate nd x nd
		m_Kt = new double[nd * nd];
		m_Ct = new double[nd * nd];
		m_Mt = new double[nd * nd];
		//allocate ns
		m_l_data = new double[ns];
		m_w_data = new double[ns];
		m_z_data = new double[nd * nz * ns];
		//allocate nd x nz
		double** data_nd_nz[] = {
			&m_r, &m_dz, &m_bt, &m_dz0r, &m_dz0t, &m_ddzr, &m_ddzt, &m_z_old, &m_z_new
		};
		m_At = new double[nd * nd * nz * nz];
		for(double** ptr : data_nd_nz) *ptr = new double[nd * nz];
	}

	//state
	void harmonic::compute_jerk(double t)
	{
		//data
		const double w = m_w_new;
		memset(m_j, 0, m_size * sizeof(double));
		//jerk
		for(uint32_t i = 1; i <= m_harmonics; i++)
		{
			const double ci = cos(i * w * t);
			const double si = sin(i * w * t);
			const double* ai = m_z_new + (2 * i - 1) * m_size;
			const double* bi = m_z_new + (2 * i + 0) * m_size;
			for(uint32_t j = 0; j < m_size; j++)
			{
				m_j[j] += i * i * i * w * w * w * si * ai[j];
				m_j[j] -= i * i * i * w * w * w * ci * bi[j];
			}
		}
	}
	void harmonic::compute_state(double t)
	{
		//data
		const double w = m_w_new;
		memcpy(m_d, m_z_new, m_size * sizeof(double));
		//state
		for(uint32_t i = 1; i <= m_harmonics; i++)
		{
			const double ci = cos(i * w * t);
			const double si = sin(i * w * t);
			const double* ai = m_z_new + (2 * i - 1) * m_size;
			const double* bi = m_z_new + (2 * i + 0) * m_size;
			for(uint32_t j = 0; j < m_size; j++)
			{
				m_d[j] += ci * ai[j];
				m_d[j] += si * bi[j];
			}
		}
	}
	void harmonic::compute_residue(void)
	{
		//data
		vector fr(m_fr, m_size);
		const vector a(m_a, m_size);
		const vector fi(m_fi, m_size);
		const vector fe(m_fe, m_size);
		const matrix Mt(m_Mt, m_size, m_size);
		//residue
		fr = m_l_new * fe - fi - Mt * a;
	}
	void harmonic::compute_velocity(double t)
	{
		//data
		const double w = m_w_new;
		memset(m_v, 0, m_size * sizeof(double));
		//velocity
		for(uint32_t i = 1; i <= m_harmonics; i++)
		{
			const double ci = cos(i * w * t);
			const double si = sin(i * w * t);
			const double* ai = m_z_new + (2 * i - 1) * m_size;
			const double* bi = m_z_new + (2 * i + 0) * m_size;
			for(uint32_t j = 0; j < m_size; j++)
			{
				m_v[j] -= i * w * si * ai[j];
				m_v[j] += i * w * ci * bi[j];
			}
		}
	}
	void harmonic::compute_acceleration(double t)
	{
		//data
		const double w = m_w_new;
		memset(m_a, 0, m_size * sizeof(double));
		//acceleration
		for(uint32_t i = 1; i <= m_harmonics; i++)
		{
			const double ci = cos(i * w * t);
			const double si = sin(i * w * t);
			const double* ai = m_z_new + (2 * i - 1) * m_size;
			const double* bi = m_z_new + (2 * i + 0) * m_size;
			for(uint32_t j = 0; j < m_size; j++)
			{
				m_a[j] -= i * i * w * w * ci * ai[j];
				m_a[j] -= i * i * w * w * si * bi[j];
			}
		}
	}
	void harmonic::compute_frequency_gradients(double t)
	{
		//data
		const double w = m_w_new;
		const uint32_t ns = m_size;
		//gradients
		vector(m_ddw, ns) = t / w * vector(m_v, ns);
		vector(m_dvw, ns) = t / w * vector(m_a, ns) + vector(m_v, ns) / w;
		vector(m_daw, ns) = t / w * vector(m_j, ns) + 2 * vector(m_a, ns) / w;
	}

	//system
	void harmonic::compute_system_residue(void)
	{
		//data
		const double w = m_w_new;
		const double T = 2 * M_PI / w;
		//size
		const uint32_t ns = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//compute
		vector fr(m_fr, ns);
		vector(m_r, ns * nz).zeros();
		for(uint32_t k = 0; k < m_quadrature_order; k++)
		{
			//quadrature
			const double sk = m_sq[k];
			const double wk = m_wq[k];
			const double tk = (1 + sk) * T / 2;
			//state
			compute_state(tk);
			compute_velocity(tk);
			compute_acceleration(tk);
			m_inertia(m_Mt, m_d, m_args);
			m_internal_force(m_fi, m_d, m_v, m_args);
			m_external_force(m_fe, tk, w, m_d, m_args);
			//residue
			compute_residue();
			vector(m_r, m_size) += wk * fr;
			for(uint32_t i = 1; i <= m_harmonics; i++)
			{
				const double ci = cos(i * w * tk);
				const double si = sin(i * w * tk);
				const uint32_t p1 = (2 * i - 1) * ns;
				const uint32_t p2 = (2 * i + 0) * ns;
				vector(m_r + p1, m_size) += wk * ci * fr;
				vector(m_r + p2, m_size) += wk * si * fr;
			}
		}
	}
	void harmonic::compute_system_tangent_p(void)
	{
		if(m_parameter == harmonic_parameter::load) compute_system_tangent_l();
		if(m_parameter == harmonic_parameter::frequency) compute_system_tangent_w();
	}
	void harmonic::compute_system_tangent_l(void)
	{
		//data
		const double w = m_w_new;
		const double T = 2 * M_PI / w;
		//size
		const uint32_t ns = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//compute
		vector fe(m_fe, ns);
		vector(m_bt, ns * nz).zeros();
		for(uint32_t k = 0; k < m_quadrature_order; k++)
		{
			//quadrature
			const double sk = m_sq[k];
			const double wk = m_wq[k];
			const double tk = (1 + sk) * T / 2;
			//state
			compute_state(tk);
			compute_velocity(tk);
			compute_acceleration(tk);
			m_external_force(m_fe, tk, m_w_new, m_d, m_args);
			//tangent
			vector(m_bt, ns) += wk * fe;
			for(uint32_t i = 1; i <= m_harmonics; i++)
			{
				const double ci = cos(i * w * tk);
				const double si = sin(i * w * tk);
				const uint32_t p1 = (2 * i - 1) * ns;
				const uint32_t p2 = (2 * i + 0) * ns;
				vector(m_bt + p1, ns) += wk * ci * fe;
				vector(m_bt + p2, ns) += wk * si * fe;
			}
		}
	}
	void harmonic::compute_system_tangent_w(void)
	{
		//data
		const double w = m_w_new;
		const double T = 2 * M_PI / w;
		//size
		const uint32_t ns = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//tangent
		vector fr(m_fr, ns);
		vector dfrw(m_dfrw, ns);
		const vector ddw(m_ddw, ns);
		const vector dvw(m_dvw, ns);
		const vector daw(m_daw, ns);
		const vector dfew(m_dfew, ns);
		const matrix Kt(m_Kt, ns, ns);
		const matrix Ct(m_Ct, ns, ns);
		const matrix Mt(m_Mt, ns, ns);
		//compute
		vector(m_bt, ns * nz).zeros();
		for(uint32_t k = 0; k < m_quadrature_order; k++)
		{
			//quadrature
			const double sk = m_sq[k];
			const double wk = m_wq[k];
			const double tk = (1 + sk) * T / 2;
			//state
			compute_jerk(tk);
			compute_state(tk);
			compute_velocity(tk);
			compute_acceleration(tk);
			m_inertia(m_Mt, m_d, m_args);
			m_damping(m_Ct, m_d, m_v, m_args);
			m_external_force(m_fe, tk, m_w_new, m_d, m_args);
			m_external_force_gradient(m_dfew, tk, m_w_new, m_d, m_args);
			m_stiffness(m_Kt, tk, m_w_new, m_l_new, m_d, m_v, m_a, m_args);
			//residue
			compute_residue();
			compute_frequency_gradients(tk);
			dfrw = m_l_new * dfew - Kt * ddw - Ct * dvw - Mt * daw;
			//tangent
			vector(m_bt, ns) += wk * dfrw;
			for(uint32_t i = 1; i <= m_harmonics; i++)
			{
				const double ci = cos(i * w * tk);
				const double si = sin(i * w * tk);
				const uint32_t p1 = (2 * i - 1) * ns;
				const uint32_t p2 = (2 * i + 0) * ns;
				vector(m_bt + p1, ns) += wk * (ci * dfrw - i * tk * si * fr);
				vector(m_bt + p2, ns) += wk * (si * dfrw + i * tk * ci * fr);
			}
		}
	}
	void harmonic::compute_system_tangent_z(void)
	{
		//data
		const double w = m_w_new;
		const double T = 2 * M_PI / w;
		//size
		const uint32_t ns = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//tangent
		matrix Kj(ns, ns);
		const matrix Kt(m_Kt, ns, ns);
		const matrix Ct(m_Ct, ns, ns);
		const matrix Mt(m_Mt, ns, ns);
		matrix At(m_At, ns * nz, ns * nz);
		//compute
		At.zeros();
		for(uint32_t k = 0; k < m_quadrature_order; k++)
		{
			//quadrature
			const double sk = m_sq[k];
			const double wk = m_wq[k];
			const double tk = (1 + sk) * T / 2;
			//residue
			compute_state(tk);
			compute_velocity(tk);
			compute_acceleration(tk);
			m_inertia(m_Mt, m_d, m_args);
			m_damping(m_Ct, m_d, m_v, m_args);
			m_stiffness(m_Kt, tk, 0, 0, m_d, m_v, m_a, m_args);
			//tangent
			At.span(0, 0, ns, ns) += wk * Kt;
			for(uint32_t i = 1; i <= m_harmonics; i++)
			{
				const double ci = cos(i * w * tk);
				const double si = sin(i * w * tk);
				const uint32_t p1 = (2 * i - 1) * ns;
				const uint32_t p2 = (2 * i + 0) * ns;
				At.span(p1, 0, ns, ns) += wk * ci * Kt;
				At.span(p2, 0, ns, ns) += wk * si * Kt;
				At.span(0, p1, ns, ns) += wk * (ci * (Kt - i * i * w * w * Mt) - i * w * si * Ct);
				At.span(0, p2, ns, ns) += wk * (si * (Kt - i * i * w * w * Mt) + i * w * ci * Ct);
				for(uint32_t j = 1; j <= m_harmonics; j++)
				{
					const double cj = cos(j * w * tk);
					const double sj = sin(j * w * tk);
					const uint32_t q1 = (2 * j - 1) * ns;
					const uint32_t q2 = (2 * j + 0) * ns;
					At.span(p1, q1, ns, ns) += wk * ci * (cj * (Kt - j * j * w * w * Mt) - j * w * sj * Ct);
					At.span(p1, q2, ns, ns) += wk * ci * (sj * (Kt - j * j * w * w * Mt) + j * w * cj * Ct);
					At.span(p2, q1, ns, ns) += wk * si * (cj * (Kt - j * j * w * w * Mt) - j * w * sj * Ct);
					At.span(p2, q2, ns, ns) += wk * si * (sj * (Kt - j * j * w * w * Mt) + j * w * cj * Ct);
				}
			}
		}
	}

	//solver
	void harmonic::record(void)
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//record
		m_l_data[m_step] = m_l_new;
		m_w_data[m_step] = m_w_new;
		memcpy(m_z_data + m_step * nd * nz, m_z_new, nd * nz * sizeof(double));
	}
	void harmonic::update(void)
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//restore
		m_l_old = m_l_new;
		m_w_old = m_w_new;
		memcpy(m_z_old, m_z_new, nd * nz * sizeof(double));
	}
	void harmonic::restore(void)
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//restore
		m_l_new = m_l_old;
		m_w_new = m_w_old;
		memcpy(m_z_new, m_z_old, nd * nz * sizeof(double));
	}

	void harmonic::increment_state(void)
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		double& p_new = m_parameter == harmonic_parameter::load ? m_l_new : m_w_new;
		const double& p_old = m_parameter == harmonic_parameter::load ? m_l_old : m_w_old;
		//increment
		p_new = p_old + m_dp;
		for(uint32_t i = 0; i < nd * nz; i++) m_z_new[i] = m_z_old[i] + m_dz[i];
	}
	void harmonic::compute_predictor(void)
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//data
		matrix At(m_At, nd * nz, nd * nz);
		vector r(m_r, nd * nz), bt(m_bt, nd * nz);
		vector dz(m_dz, nd * nz), dz0r(m_dz0r, nd * nz), dz0t(m_dz0t, nd * nz);
		//state
		compute_system_residue();
		compute_system_tangent_p();
		compute_system_tangent_z();
		//parameter
		At.solve(dz0r, r);
		At.solve(dz0t, bt);
		m_dp = m_step == 0 ? m_dpg : compute_parameter_predictor();
		//state
		dz = dz0r + m_dp * dz0t;
		//increment
		increment_state();
	}
	void harmonic::compute_corrector(void)
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//data
		matrix At(m_At, nd * nz, nd * nz);
		vector r(m_r, nd * nz), bt(m_bt, nd * nz);
		vector dz(m_dz, nd * nz), ddzr(m_ddzr, nd * nz), ddzt(m_ddzt, nd * nz);
		//state
		compute_system_residue();
		compute_system_tangent_p();
		compute_system_tangent_z();
		//corrector
		At.solve(ddzr, r);
		At.solve(ddzt, bt);
		m_ddp = compute_parameter_corrector();
		//state
		m_dp += m_ddp;
		dz += ddzr + m_ddp * ddzt;
		//increment
		increment_state();
	}

	double harmonic::compute_parameter_predictor(void) const
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//data
		const vector dz(m_dz, nd * nz);
		const vector dz0r(m_dz0r, nd * nz);
		const vector dz0t(m_dz0t, nd * nz);
		//predictor
		const double a = dz0t.inner(dz0t);
		const double b = dz0t.inner(dz0r);
		const double s = math::sign(dz0t.inner(dz));
		const double c = dz0r.inner(dz0r) - dz.inner(dz);
		//return
		return -b / a + s * sqrt(b * b - a * c) / a;
	}
	double harmonic::compute_parameter_corrector(void) const
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//data
		const vector dz(m_dz, nd * nz);
		const vector ddzr(m_ddzr, nd * nz);
		const vector ddzt(m_ddzt, nd * nz);
		//corrector
		const double a = ddzt.inner(ddzt);
		const double b = ddzt.inner(ddzr + dz);
		const double c = ddzr.inner(ddzr + 2 * dz);
		const double s = math::sign(ddzt.inner(dz));
		//return
		return -b / a + s * sqrt(b * b - a * c) / a;
	}
}