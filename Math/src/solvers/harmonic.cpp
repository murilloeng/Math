//math
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/solvers/harmonic.hpp"

namespace math
{
	//constructors
	harmonic::harmonic(void) : 
		m_args(nullptr), m_size(0), m_harmonics(0), m_iteration_max(0), m_quadrature_order(0),
		m_frequency(0.0), m_tolerance(0.0), 
		m_external_force(nullptr), m_internal_force(nullptr), m_inertia(nullptr), m_damping(nullptr), m_stiffness(nullptr),
		m_quadrature(nullptr), m_d(nullptr), m_v(nullptr), m_a(nullptr), m_r(nullptr),
		m_fi(nullptr), m_fe(nullptr), m_Kt(nullptr), m_Ct(nullptr), m_Mt(nullptr), m_At(nullptr)
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
		//setup
		vector dz(m_size * (2 * m_harmonics + 1));
		vector(m_z, m_size * (2 * m_harmonics + 1)).zeros();
		matrix At(m_At, m_size * (2 * m_harmonics + 1), m_size * (2 * m_harmonics + 1));
		//solve
		for(uint32_t iteration = 0; iteration < m_iteration_max; iteration++)
		{
			//check
			compute_residue();
			if(vector(m_r, m_size * (2 * m_harmonics + 1)).norm() < m_tolerance)
			{
				return true;
			}
			//update
			compute_tangent();
			At.solve(dz, vector(m_r, m_size * (2 * m_harmonics + 1)));
			vector(m_z, m_size * (2 * m_harmonics + 1)) += dz;
		}
		//return
		return false;
	}
	void harmonic::setup(void)
	{
		cleanup();
		m_d = new double[m_size];
		m_v = new double[m_size];
		m_a = new double[m_size];
		m_fi = new double[m_size];
		m_fe = new double[m_size];
		m_fr = new double[m_size];
		m_Kt = new double[m_size * m_size];
		m_Ct = new double[m_size * m_size];
		m_Mt = new double[m_size * m_size];
		m_z = new double[m_size * (2 * m_harmonics + 1)];
		m_r = new double[m_size * (2 * m_harmonics + 1)];
		m_At = new double[m_size * m_size * (2 * m_harmonics + 1) * (2 * m_harmonics + 1)];
		m_quadrature = new quadrature::Quadrature(quadrature::rule::legendre, m_quadrature_order);

	}
	void harmonic::cleanup(void)
	{
		delete[] m_d;
		delete[] m_v;
		delete[] m_a;
		delete[] m_r;
		delete[] m_fi;
		delete[] m_fe;
		delete[] m_fr;
		delete[] m_Kt;
		delete[] m_Ct;
		delete[] m_Mt;
		delete[] m_At;
		delete m_quadrature;
	}

	//solver
	void harmonic::compute_residue(void)
	{
		//data
		const double w = m_frequency;
		const double T = 2 * M_PI / w;
		memset(m_r, 0, m_size * (2 * m_harmonics + 1) * sizeof(double));
		//compute
		for(uint32_t k = 0; k < m_quadrature_order; k++)
		{
			//quadrature
			const double sk = m_quadrature->point(k);
			const double wk = m_quadrature->weight(k);
			const double tk = (1 + sk) * T / 2;
			//residue
			compute_residue(tk);
			vector(m_r, m_size) += wk * vector(m_fr, m_size);
			for(uint32_t i = 1; i <= m_harmonics; i++)
			{
				const double ci = cos(i * w * tk);
				const double si = sin(i * w * tk);
				vector(m_r + (2 * i - 1) * m_size, m_size) += wk * ci * vector(m_fr, m_size);
				vector(m_r + (2 * i + 0) * m_size, m_size) += wk * si * vector(m_fr, m_size);
			}
		}
	}
	void harmonic::compute_tangent(void)
	{
		//data
		const double w = m_frequency;
		const double T = 2 * M_PI / w;
		//tangent
		const uint32_t ns = m_size;
		const uint32_t nh = m_harmonics;
		const matrix Kt(m_Kt, m_size, m_size);
		const matrix Ct(m_Ct, m_size, m_size);
		const matrix Mt(m_Mt, m_size, m_size);
		matrix At(m_At, ns * (2 * nh + 1), ns * (2 * nh + 1));
		//compute
		At.zeros();
		for(uint32_t k = 0; k < m_quadrature_order; k++)
		{
			//quadrature
			const double sk = m_quadrature->point(k);
			const double wk = m_quadrature->weight(k);
			const double tk = (1 + sk) * T / 2;
			//residue
			compute_residue(tk);
			m_damping(m_Ct, m_d, m_v, m_args);
			m_stiffness(m_Kt, tk, m_d, m_v, m_args);
			//tangent
			At.span(0, 0, ns, ns) += wk * Kt;
			for(uint32_t i = 1; i <= nh; i++)
			{
				const double ci = cos(i * w * tk);
				const double si = sin(i * w * tk);
				At.span((2 * i - 1) * ns, 0, ns, ns) += wk * ci * Kt;
				At.span((2 * i + 0) * ns, 0, ns, ns) += wk * si * Kt;
				At.span(0, (2 * i - 1) * ns, ns, ns) += wk * (ci * (Kt - i * i * w * w * Mt) - i * w * si * Ct);
				At.span(0, (2 * i + 0) * ns, ns, ns) += wk * (si * (Kt - i * i * w * w * Mt) + i * w * ci * Ct);
				for(uint32_t j = 1; j <= nh; j++)
				{
					const double cj = cos(j * w * tk);
					const double sj = sin(j * w * tk);
					At.span((2 * i - 1) * ns, (2 * j - 1) * ns, ns, ns) += wk * ci * (cj * (Kt - j * j * w * w * Mt) - j * w * sj * Ct);
					At.span((2 * i - 1) * ns, (2 * j + 0) * ns, ns, ns) += wk * ci * (sj * (Kt - j * j * w * w * Mt) + j * w * cj * Ct);
					At.span((2 * i + 0) * ns, (2 * j - 1) * ns, ns, ns) += wk * si * (cj * (Kt - j * j * w * w * Mt) - j * w * sj * Ct);
					At.span((2 * i + 0) * ns, (2 * j + 0) * ns, ns, ns) += wk * si * (sj * (Kt - j * j * w * w * Mt) + j * w * cj * Ct);
				}
			}
		}
	}

	void harmonic::compute_state(double t)
	{
		const double w = m_frequency;
		memcpy(m_d, m_z, m_size * sizeof(double));
		for(uint32_t i = 1; i <= m_harmonics; i++)
		{
			const double c = cos(i * w * t);
			const double s = sin(i * w * t);
			for(uint32_t j = 0; j < m_size; j++)
			{
				m_d[j] += c * m_z[(2 * i - 1) * m_size + j];
				m_d[j] += s * m_z[(2 * i + 0) * m_size + j];
			}
		}
	}
	void harmonic::compute_residue(double t)
	{
		//data
		vector fr(m_fr, m_size);
		const vector a(m_a, m_size);
		const vector fi(m_fi, m_size);
		const vector fe(m_fe, m_size);
		const matrix Mt(m_Mt, m_size, m_size);
		//state
		compute_state(t);
		compute_velocity(t);
		compute_acceleration(t);
		m_inertia(m_Mt, m_d, m_args);
		m_external_force(m_fe, t, m_d, m_args);
		m_internal_force(m_fi, m_d, m_v, m_args);
		//compute
		fr = fe - fi - Mt * a;
	}
	void harmonic::compute_velocity(double t)
	{
		const double w = m_frequency;
		memset(m_v, 0, m_size * sizeof(double));
		for(uint32_t i = 1; i <= m_harmonics; i++)
		{
			const double c = cos(i * w * t);
			const double s = sin(i * w * t);
			for(uint32_t j = 0; j < m_size; j++)
			{
				m_v[j] -= i * w * s * m_z[(2 * i - 1) * m_size + j];
				m_v[j] += i * w * c * m_z[(2 * i + 0) * m_size + j];
			}
		}
	}
	void harmonic::compute_acceleration(double t)
	{
		const double w = m_frequency;
		memset(m_a, 0, m_size * sizeof(double));
		for(uint32_t i = 1; i <= m_harmonics; i++)
		{
			const double c = cos(i * w * t);
			const double s = sin(i * w * t);
			for(uint32_t j = 0; j < m_size; j++)
			{
				m_a[j] -= i * i * w * w * c * m_z[(2 * i - 1) * m_size + j];
				m_a[j] -= i * i * w * w * s * m_z[(2 * i + 0) * m_size + j];
			}
		}
	}
}