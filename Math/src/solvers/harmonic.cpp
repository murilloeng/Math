//std
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <cstring>
#include <malloc.h>

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
		m_args(nullptr), m_size(0), m_step_max(0), m_harmonics(0), 
		m_attempt_max(0), m_iteration_max(0), m_quadrature_order(0), 
		m_dpg(0), m_tolerance(0), 
		m_l_0(0), m_l_min(-DBL_MAX), m_l_max(+DBL_MAX), 
		m_w_0(0), m_w_min(-DBL_MAX), m_w_max(+DBL_MAX),
		m_strategy(harmonic_strategy::uniform_increment), m_control(harmonic_control::load),
		m_internal_force(nullptr), m_external_force(nullptr), 
		m_inertia(nullptr), m_damping(nullptr), m_stiffness(nullptr),
		m_sq(nullptr), m_wq(nullptr), 
		m_d(nullptr), m_v(nullptr), m_a(nullptr), m_dvw(nullptr), m_daw(nullptr),
		m_l_data(nullptr), m_w_data(nullptr), m_z_old(nullptr), m_z_new(nullptr), m_z_data(nullptr),
		m_dp(0), m_ddp(0),
		m_dz(nullptr), m_dz0r(nullptr), m_dz0t(nullptr), m_ddzr(nullptr), m_ddzt(nullptr),
		m_r(nullptr), m_fi(nullptr), m_fe(nullptr), m_fr(nullptr),
		m_Kt(nullptr), m_Ct(nullptr), m_Mt(nullptr), m_At(nullptr), m_bt(nullptr), m_dfrw(nullptr),
		m_stability_data(nullptr), m_y(nullptr), m_y1(nullptr), m_y2(nullptr), m_y3(nullptr), m_y4(nullptr)
	{
		return;
	}

	//destructor
	harmonic::~harmonic(void)
	{
		cleanup();
	}

	//solve
	void harmonic::save(void)
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//open
		FILE* file = fopen("harmonic.txt", "w");
		//save
		for(uint32_t i = 0; i <= m_step; i++)
		{
			fprintf(file, "%+.6e ", m_l_data[i]);
			fprintf(file, "%+.6e ", m_w_data[i]);
			for(uint32_t j = 0; j < nd * nz; j++)
			{
				fprintf(file, "%+.6e ", m_z_data[i * nd * nz + j]);
			}
			fprintf(file, "\n");
		}
		//close
		fclose(file);
	}
	bool harmonic::solve(void)
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//solve
		setup();
		record(0);
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
				update();
				record(m_step + 1);
				printf("step: %04d iteration: %02d\n", m_step, m_iteration);
			}
			if(stop()) break;
		}
		//return
		return true;
	}

	//solver
	bool harmonic::stop(void)
	{
		if(m_control == harmonic_control::load && m_l_new < m_l_min)
		{
			printf("Min load reached!\n");
			return true;
		}
		if(m_control == harmonic_control::load && m_l_new > m_l_max)
		{
			printf("Max load reached!\n");
			return true;
		}
		if(m_control == harmonic_control::frequency && m_w_new < m_w_min)
		{
			printf("Min frequency reached!\n");
			return true;
		}
		if(m_control == harmonic_control::frequency && m_w_new > m_w_max)
		{
			printf("Max frequency reached!\n");
			return true;
		}
		return false;
	}
	void harmonic::setup(void)
	{
		cleanup();
		allocate();
		initialize();
	}
	void harmonic::cleanup(void)
	{
		double* data[] = {
			m_sq, m_wq,
			m_d, m_v, m_a, m_dvw, m_daw,
			m_l_data, m_w_data, m_z_old, m_z_new, m_z_data,
			m_dz, m_dz0r, m_dz0t, m_ddzr, m_ddzt,
			m_r, m_fi, m_fe, m_fr, m_Kt, m_Ct, m_Mt, m_At, m_bt, m_dfrw,
			m_y, m_y1, m_y2, m_y3, m_y4
		};
		delete[] m_stability_data;
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
			&m_d, &m_v, &m_a, &m_fi, &m_fe, &m_fr, &m_dvw, &m_daw, &m_dfrw
		};
		for(double** ptr : data_nd) *ptr = new double[nd];
		//allocate nq
		m_sq = new double[nq];
		m_wq = new double[nq];
		//allocate ns
		m_l_data = new double[ns + 1];
		m_w_data = new double[ns + 1];
		m_stability_data = new bool[ns + 1];
		m_z_data = new double[nd * nz * (ns + 1)];
		//allocate nd x nd
		double** data_nd_nd[] = {
			&m_y, &m_y1, &m_y2, &m_y3, &m_y4, &m_Kt, &m_Ct, &m_Mt
		};
		for(double** ptr : data_nd_nd) *ptr = new double[nd * nd];
		//allocate nd x nz
		double** data_nd_nz[] = {
			&m_r, &m_dz, &m_bt, &m_dz0r, &m_dz0t, &m_ddzr, &m_ddzt, &m_z_old, &m_z_new
		};
		m_At = new double[nd * nd * nz * nz];
		for(double** ptr : data_nd_nz) *ptr = new double[nd * nz];
	}
	void harmonic::initialize(void)
	{
		//data
		bool equilibrium;
		const uint32_t nd = m_size;
		const uint32_t nq = m_quadrature_order;
		const uint32_t nz = 2 * m_harmonics + 1;
		//quadrature
		legendre_dr_compute(m_quadrature_order, m_sq, m_wq);
		//state
		m_l_old = m_l_new = m_l_0;
		m_w_old = m_w_new = m_w_0;
		matrix At(m_At, nd * nz, nd * nz);
		vector r(m_r, nd * nz), dz(m_dz, nd * nz);
		memset(m_z_new, 0, nd * nz * sizeof(double));
		for(m_iteration = 0; m_iteration < m_iteration_max; m_iteration++)
		{
			//residue
			compute_system_residue();
			compute_system_tangent_z();
			equilibrium = r.norm() < m_tolerance;
			//increment
			At.solve(dz, r);
			if(equilibrium) break;
			vector(m_z_new, nd * nz) += dz;
		}
		if(!equilibrium)
		{
			printf("Harmonic solver failed to converge at initialization\n");
			return;
		}
		memcpy(m_z_old, m_z_new, nd * nz * sizeof(double));
	}

	//state
	void harmonic::compute_state(double t)
	{
		//data
		const double w = m_w_new;
		const uint32_t nd = m_size;
		memcpy(m_d, m_z_new, nd * sizeof(double));
		//state
		for(uint32_t i = 1; i <= m_harmonics; i++)
		{
			const double ci = cos(i * w * t);
			const double si = sin(i * w * t);
			const uint32_t p1 = (2 * i - 1) * nd;
			const uint32_t p2 = (2 * i + 0) * nd;
			for(uint32_t j = 0; j < nd; j++)
			{
				m_d[j] += ci * m_z_new[p1 + j];
				m_d[j] += si * m_z_new[p2 + j];
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
		const uint32_t nd = m_size;
		memset(m_v, 0, nd * sizeof(double));
		//velocity
		for(uint32_t i = 1; i <= m_harmonics; i++)
		{
			const double ci = cos(i * w * t);
			const double si = sin(i * w * t);
			const uint32_t p1 = (2 * i - 1) * nd;
			const uint32_t p2 = (2 * i + 0) * nd;
			for(uint32_t j = 0; j < nd; j++)
			{
				m_v[j] -= i * w * si * m_z_new[p1 + j];
				m_v[j] += i * w * ci * m_z_new[p2 + j];
			}
		}
	}
	void harmonic::compute_acceleration(double t)
	{
		//data
		const double w = m_w_new;
		const uint32_t nd = m_size;
		memset(m_a, 0, nd * sizeof(double));
		//acceleration
		for(uint32_t i = 1; i <= m_harmonics; i++)
		{
			const double ci = cos(i * w * t);
			const double si = sin(i * w * t);
			const uint32_t p1 = (2 * i - 1) * nd;
			const uint32_t p2 = (2 * i + 0) * nd;
			for(uint32_t j = 0; j < nd; j++)
			{
				m_a[j] -= i * i * w * w * ci * m_z_new[p1 + j];
				m_a[j] -= i * i * w * w * si * m_z_new[p2 + j];
			}
		}
	}
	void harmonic::compute_frequency_gradients(double t)
	{
		//data
		const double w = m_w_new;
		const uint32_t nd = m_size;
		memset(m_dvw, 0, nd * sizeof(double));
		memset(m_daw, 0, nd * sizeof(double));
		//gradients
		for(uint32_t i = 1; i <= m_harmonics; i++)
		{
			const double ci = cos(i * w * t);
			const double si = sin(i * w * t);
			const uint32_t p1 = (2 * i - 1) * nd;
			const uint32_t p2 = (2 * i + 0) * nd;
			for(uint32_t j = 0; j < nd; j++)
			{
				m_dvw[j] -= i * si * m_z_new[p1 + j];
				m_dvw[j] += i * ci * m_z_new[p2 + j];
				m_daw[j] -= 2 * i * i * w * ci * m_z_new[p1 + j];
				m_daw[j] -= 2 * i * i * w * si * m_z_new[p2 + j];
			}
		}
	}

	//system
	void harmonic::compute_system_residue(void)
	{
		//data
		const double w = m_w_new;
		const double T = 2 * M_PI / w;
		//size
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//compute
		vector fr(m_fr, nd);
		vector(m_r, nd * nz).zeros();
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
			m_external_force(m_fe, tk, m_w_new, m_d, m_args);
			//residue
			compute_residue();
			vector(m_r, m_size) += wk * fr;
			for(uint32_t i = 1; i <= m_harmonics; i++)
			{
				const double ci = cos(i * w * tk);
				const double si = sin(i * w * tk);
				const uint32_t p1 = (2 * i - 1) * nd;
				const uint32_t p2 = (2 * i + 0) * nd;
				vector(m_r + p1, m_size) += wk * ci * fr;
				vector(m_r + p2, m_size) += wk * si * fr;
			}
		}
	}
	void harmonic::compute_system_tangent_p(void)
	{
		if(m_control == harmonic_control::load) compute_system_tangent_l();
		if(m_control == harmonic_control::frequency) compute_system_tangent_w();
	}
	void harmonic::compute_system_tangent_l(void)
	{
		//data
		const double w = m_w_new;
		const double T = 2 * M_PI / w;
		//size
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//compute
		vector fe(m_fe, nd);
		vector(m_bt, nd * nz).zeros();
		for(uint32_t k = 0; k < m_quadrature_order; k++)
		{
			//quadrature
			const double sk = m_sq[k];
			const double wk = m_wq[k];
			const double tk = (1 + sk) * T / 2;
			//state
			compute_state(tk);
			m_external_force(m_fe, tk, w, m_d, m_args);
			//tangent
			vector(m_bt, nd) += wk * fe;
			for(uint32_t i = 1; i <= m_harmonics; i++)
			{
				const double ci = cos(i * w * tk);
				const double si = sin(i * w * tk);
				const uint32_t p1 = (2 * i - 1) * nd;
				const uint32_t p2 = (2 * i + 0) * nd;
				vector(m_bt + p1, nd) += wk * ci * fe;
				vector(m_bt + p2, nd) += wk * si * fe;
			}
		}
	}
	void harmonic::compute_system_tangent_w(void)
	{
		//data
		const double w = m_w_new;
		const double T = 2 * M_PI / w;
		//size
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//tangent
		vector dfrw(m_dfrw, nd);
		const vector dvw(m_dvw, nd);
		const vector daw(m_daw, nd);
		const matrix Kt(m_Kt, nd, nd);
		const matrix Ct(m_Ct, nd, nd);
		const matrix Mt(m_Mt, nd, nd);
		//compute
		vector(m_bt, nd * nz).zeros();
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
			m_damping(m_Ct, m_d, m_v, m_args);
			//gradient
			compute_frequency_gradients(tk);
			dfrw = -Ct * dvw - Mt * daw;
			//tangent
			vector(m_bt, nd) += wk * dfrw;
			for(uint32_t i = 1; i <= m_harmonics; i++)
			{
				const double ci = cos(i * w * tk);
				const double si = sin(i * w * tk);
				const uint32_t p1 = (2 * i - 1) * nd;
				const uint32_t p2 = (2 * i + 0) * nd;
				vector(m_bt + p1, nd) += wk * ci * dfrw;
				vector(m_bt + p2, nd) += wk * si * dfrw;
			}
		}
	}
	void harmonic::compute_system_tangent_z(void)
	{
		//data
		const double w = m_w_new;
		const double T = 2 * M_PI / w;
		//size
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//tangent
		const matrix Kt(m_Kt, nd, nd);
		const matrix Ct(m_Ct, nd, nd);
		const matrix Mt(m_Mt, nd, nd);
		matrix At(m_At, nd * nz, nd * nz);
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
			m_stiffness(m_Kt, tk, w, m_l_new, m_d, m_v, m_a, m_args);
			//tangent
			At.span(0, 0, nd, nd) += wk * Kt;
			for(uint32_t i = 1; i <= m_harmonics; i++)
			{
				const double ci = cos(i * w * tk);
				const double si = sin(i * w * tk);
				const uint32_t p1 = (2 * i - 1) * nd;
				const uint32_t p2 = (2 * i + 0) * nd;
				At.span(p1, 0, nd, nd) += wk * ci * Kt;
				At.span(p2, 0, nd, nd) += wk * si * Kt;
				At.span(0, p1, nd, nd) += wk * (ci * (Kt - i * i * w * w * Mt) - i * w * si * Ct);
				At.span(0, p2, nd, nd) += wk * (si * (Kt - i * i * w * w * Mt) + i * w * ci * Ct);
				for(uint32_t j = 1; j <= m_harmonics; j++)
				{
					const double cj = cos(j * w * tk);
					const double sj = sin(j * w * tk);
					const uint32_t q1 = (2 * j - 1) * nd;
					const uint32_t q2 = (2 * j + 0) * nd;
					At.span(p1, q1, nd, nd) += wk * ci * (cj * (Kt - j * j * w * w * Mt) - j * w * sj * Ct);
					At.span(p1, q2, nd, nd) += wk * ci * (sj * (Kt - j * j * w * w * Mt) + j * w * cj * Ct);
					At.span(p2, q1, nd, nd) += wk * si * (cj * (Kt - j * j * w * w * Mt) - j * w * sj * Ct);
					At.span(p2, q2, nd, nd) += wk * si * (sj * (Kt - j * j * w * w * Mt) + j * w * cj * Ct);
				}
			}
		}
	}

	//solver
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
	void harmonic::record(uint32_t step)
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//record
		m_l_data[step] = m_l_new;
		m_w_data[step] = m_w_new;
		memcpy(m_z_data + step * nd * nz, m_z_new, nd * nz * sizeof(double));
	}

	void harmonic::increment_state(void)
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		double& p_new = m_control == harmonic_control::load ? m_l_new : m_w_new;
		const double& p_old = m_control == harmonic_control::load ? m_l_old : m_w_old;
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
		compute_parameter_predictor();
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
		compute_parameter_corrector();
		//state
		m_dp += m_ddp;
		dz += ddzr + m_ddp * ddzt;
		//increment
		increment_state();
	}

	void harmonic::compute_parameter_predictor(void)
	{
		if(m_step == 0)
		{
			m_dp = m_dpg;
		}
		else
		{
			if(m_strategy == harmonic_strategy::minimal_norm)
				predictor_minimal_norm();
			if(m_strategy == harmonic_strategy::uniform_increment)
				predictor_uniform_increment();
			if(m_strategy == harmonic_strategy::arc_length_spherical)
				predictor_arc_length_spheric();
			if(m_strategy == harmonic_strategy::arc_length_cylindrical)
				predictor_arc_length_cylindric();
			if(isnan(m_dp)) printf("parameter predictor failed!\n");
		}
	}
	void harmonic::compute_parameter_corrector(void)
	{
		if(m_strategy == harmonic_strategy::minimal_norm)
			corrector_minimal_norm();
		if(m_strategy == harmonic_strategy::uniform_increment)
			corrector_uniform_increment();
		if(m_strategy == harmonic_strategy::arc_length_spherical)
			corrector_arc_length_spheric();
		if(m_strategy == harmonic_strategy::arc_length_cylindrical)
			corrector_arc_length_cylindric();
		if(isnan(m_ddp)) printf("parameter corrector failed!\n");
	}

	void harmonic::predictor_minimal_norm(void)
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//data
		vector dz(m_dz, nd * nz);
		vector dz0t(m_dz0t, nd * nz);
		//predictor
		m_dp = math::sign(dz.inner(dz0t)) * dz.norm() / dz0t.norm();
	}
	void harmonic::corrector_minimal_norm(void)
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//data
		vector ddzr(m_ddzr, nd * nz);
		vector ddzt(m_ddzt, nd * nz);
		//corrector
		m_ddp = -ddzt.inner(ddzr) / ddzt.inner(ddzt);
	}

	void harmonic::predictor_uniform_increment(void)
	{
		m_dp = m_dpg;
	}
	void harmonic::corrector_uniform_increment(void)
	{
		m_ddp = 0;
	}

	void harmonic::predictor_arc_length_spheric(void)
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//data
		const vector dz(m_dz, nd * nz);
		const vector dz0r(m_dz0r, nd * nz);
		const vector dz0t(m_dz0t, nd * nz);
		//predictor
		const double b = dz0t.inner(dz0r);
		const double a = dz0t.inner(dz0t) + 1;
		const double s = math::sign(dz0t.inner(dz));
		const double c = dz0r.inner(dz0r) - dz.inner(dz) - m_dp * m_dp;
		//predictor
		m_dp = -b / a + s * sqrt(b * b - a * c) / a;
	}
	void harmonic::corrector_arc_length_spheric(void)
	{
		//data
		const uint32_t nd = m_size;
		const uint32_t nz = 2 * m_harmonics + 1;
		//data
		const vector dz(m_dz, nd * nz);
		const vector ddzr(m_ddzr, nd * nz);
		const vector ddzt(m_ddzt, nd * nz);
		//corrector
		const double a = ddzt.inner(ddzt) + 1;
		const double c = ddzr.inner(ddzr + 2 * dz);
		const double s = math::sign(ddzt.inner(dz));
		const double b = ddzt.inner(ddzr + dz) + m_dp;
		//corrector
		m_ddp = -b / a + s * sqrt(b * b - a * c) / a;
	}

	void harmonic::predictor_arc_length_cylindric(void)
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
		//predictor
		m_dp = -b / a + s * sqrt(b * b - a * c) / a;
	}
	void harmonic::corrector_arc_length_cylindric(void)
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
		//corrector
		m_ddp = -b / a + s * sqrt(b * b - a * c) / a;
	}
}