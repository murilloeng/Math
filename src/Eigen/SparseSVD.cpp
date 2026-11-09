//std
#include <cmath>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Linear/Sparse.hpp"
#include "Math/inc/Eigen/SparseSVD.hpp"
#include "Math/inc/Miscellaneous/util.hpp"

extern "C"
{
	void dsaupd_(int32_t*, const char*, const uint32_t*, const char*, const uint32_t*, const double*, double*, const uint32_t*, double*, const uint32_t*, int32_t*, int32_t*, double*, double*, const uint32_t*, int32_t*);
	void dseupd_(const uint32_t*, const char*, uint32_t*, double*, double*, const uint32_t*, const double*, const char*, const uint32_t*, const char*, const uint32_t*, const double*, double*, const uint32_t*, double*, const uint32_t*, int32_t*, int32_t*, double*, double*, const uint32_t*, int32_t*);
}

//A = U * S * V^T is computed from the largest eigenpairs of A^T * A:
//values below sqrt(eps) * s_max are not resolved and their U columns are set to zero

namespace math
{
	namespace eigen
	{
		//constructor
		SparseSVD::SparseSVD(void) : SparseSym{0U, 0U, 0U, nullptr, nullptr}, m_rows{0U}, m_s{nullptr}, m_U{nullptr}, m_V{nullptr}, m_w{nullptr}, m_A{nullptr}
		{
			type("LA");
		}
		SparseSVD::SparseSVD(uint32_t rows, uint32_t cols, uint32_t modes, uint32_t vectors, const int32_t* rows_map, const int32_t* cols_map, const double* A, double* s, double* U, double* V) :
			SparseSym{cols, modes, vectors, rows_map, cols_map}, m_rows{rows}, m_s{s}, m_U{U}, m_V{V}, m_w{nullptr}, m_A{A}
		{
			type("LA");
		}

		//data
		uint32_t SparseSVD::rows(void) const
		{
			return m_rows;
		}
		uint32_t SparseSVD::cols(void) const
		{
			return m_order;
		}

		//compute
		bool SparseSVD::compute(void)
		{
			//data
			const uint32_t n = m_order;
			const uint32_t nev = m_modes;
			const uint32_t ncv = m_vectors;
			const uint32_t lworkl = (ncv + 8) * ncv;
			//check
			if(nev == 0 || nev > m_rows || nev >= ncv || ncv > n)
			{
				throw std::runtime_error("SparseSVD requires 0 < modes <= rows and modes < vectors <= cols!");
			}
			//data
			const char* which = m_type;
			const double tol = m_tolerance;
			int32_t info = 0, iparam[11], ipntr[11];
			double* resid = new double[n], *v = new double[n * ncv];
			double* workd = new double[3 * n], *workl = new double[lworkl];
			m_w = new double[m_rows];
			//setup
			m_ido = 0;
			iparam[0] = 1;
			iparam[6] = m_mode;
			iparam[2] = m_iteration_max;
			//Lanczos
			while(true)
			{
				//decomposition
				dsaupd_(&m_ido, "I", &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
				//compute
				if(m_ido == 99) break;
				if(m_ido == -1 || m_ido == 1) operation(workd + ipntr[1] - 1, workd + ipntr[0] - 1);
			}
			if(info == 0)
			{
				//eigenvalues
				const uint32_t rvec = bool(m_V);
				uint32_t* select = new uint32_t[ncv];
				dseupd_(&rvec, "A", select, m_s, m_V, &n, &m_shift, "I", &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
				delete[] select;
			}
			if(info == 0)
			{
				//sort
				for(uint32_t i = 0; i < nev / 2; i++)
				{
					swap(m_s[i], m_s[nev - 1 - i]);
					for(uint32_t j = 0; m_V && j < n; j++)
					{
						swap(m_V[j + n * i], m_V[j + n * (nev - 1 - i)]);
					}
				}
				//singular values
				for(uint32_t i = 0; i < nev; i++)
				{
					m_s[i] = sqrt(fmax(m_s[i], 0));
				}
				//left vectors
				const double s_min = 1.00e-8 * m_s[0];
				const math::Sparse As(m_A, m_rows_map, m_cols_map, m_rows, n);
				for(uint32_t i = 0; m_U && m_V && i < nev; i++)
				{
					As.product(m_U + m_rows * i, m_V + n * i);
					for(uint32_t j = 0; j < m_rows; j++)
					{
						m_U[j + m_rows * i] = m_s[i] > s_min ? m_U[j + m_rows * i] / m_s[i] : 0;
					}
				}
			}
			//delete
			delete[] v;
			delete[] m_w;
			delete[] resid;
			delete[] workd;
			delete[] workl;
			m_w = nullptr;
			//return
			return info == 0;
		}
		void SparseSVD::setup(uint32_t rows, uint32_t cols, uint32_t modes, uint32_t vectors, const int32_t* rows_map, const int32_t* cols_map, const double* A, double* s, double* U, double* V)
		{
			m_s = s;
			m_U = U;
			m_V = V;
			m_A = A;
			m_rows = rows;
			m_order = cols;
			m_modes = modes;
			m_vectors = vectors;
			m_rows_map = rows_map;
			m_cols_map = cols_map;
		}

		//operation
		void SparseSVD::operation(double* y, double* x) const
		{
			//w = A * x
			math::Sparse(m_A, m_rows_map, m_cols_map, m_rows, m_order).product(m_w, x);
			//y = A^T * w
			for(uint32_t i = 0; i < m_order; i++)
			{
				y[i] = 0;
				for(int32_t j = m_cols_map[i]; j < m_cols_map[i + 1]; j++)
				{
					y[i] += m_A[j] * m_w[m_rows_map[j]];
				}
			}
		}
	}
}
