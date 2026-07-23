//std
#include <cstdio>
#include <cstring>

//Math
#include "Math/inc/Linear/Sparse.hpp"
#include "Math/inc/Eigen/SparseNonGen.hpp"

extern "C"
{
	void dnaupd_(int32_t*, const char*, const uint32_t*, const char*, const uint32_t*, const double*, double*, const uint32_t*, double*, const uint32_t*, int32_t*, int32_t*, double*, double*, const uint32_t*, int32_t*);
	void dneupd_(const uint32_t*, const char*, uint32_t*, double*, double*, double*, const uint32_t*, const double*, const double*, double*, const char*, const uint32_t*, const char*, const uint32_t*, const double*, double*, const uint32_t*, double*, const uint32_t*, int32_t*, int32_t*, double*, double*, const uint32_t*, int32_t*);
}

namespace math
{
	namespace eigen
	{
		//constructor
		SparseNonGen::SparseNonGen(void) : SparseNon{0U, 0U, 0U, nullptr, nullptr}, m_U{nullptr}, m_sr{nullptr}, m_si{nullptr}, m_A{nullptr}
		{
			m_mode = 2;
		}
		SparseNonGen::SparseNonGen(uint32_t order, uint32_t modes, uint32_t vectors, const int32_t* rows_map, const int32_t* cols_map, const double* A, const double* B, double* sr, double* si, double* U) : SparseNon{order, modes, vectors, rows_map, cols_map}, m_U{U}, m_sr{sr}, m_si{si}, m_A{A}
		{
			m_mode = 2;
		}

		//destructor
		SparseNonGen::~SparseNonGen(void)
		{
			return;
		}

		//compute
		bool SparseNonGen::compute(void)
		{
			//data
			const uint32_t n = m_order;
			const uint32_t nev = m_modes;
			const uint32_t ncv = m_vectors;
			const uint32_t lworkl = (3 * ncv + 6) * ncv;
			const math::Sparse Bs(m_B, m_rows_map, m_cols_map, m_order, m_order);
			//data
			m_S = &Bs;
			m_S->decompose();
			const char* which = m_type;
			const double tol = m_tolerance;
			int32_t info = 0, iparam[11], ipntr[11];
			double *resid = new double[n], *v = new double[n * ncv];
			double *workd = new double[3 * n], *workl = new double[lworkl];
			//setup
			m_ido = 0;
			iparam[0] = 1;
			iparam[6] = m_mode;
			iparam[2] = m_iteration_max;
			//Arnoldi
			while(true)
			{
				//decomposition
				dnaupd_(&m_ido, "G", &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
				//compute
				if(m_ido == 99) break;
				if(m_ido == -1 || m_ido == 1 || m_ido == 2) operation(workd + ipntr[1] - 1, workd + ipntr[0] - 1);
			}
			if(info != 0) printf("info: %d\n", info);
			if(info != 0) return false;
			//eigenvalues
			m_S->release();
			const uint32_t rvec = bool(m_U);
			double* workev = new double[3 * n];
			uint32_t* select = new uint32_t[ncv];
			dneupd_(&rvec, "A", select, m_sr, m_si, m_U, &n, &m_shift[0], &m_shift[1], workev, "G", &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
			//delete
			delete[] v;
			delete[] resid;
			delete[] workd;
			delete[] workl;
			delete[] workev;
			delete[] select;
			//return
			return info == 0;
		}
		void SparseNonGen::setup(uint32_t order, uint32_t modes, uint32_t vectors, const int32_t* rows_map, const int32_t* cols_map, const double* A, const double* B, double* sr, double* si, double* U)
		{
			m_U = U;
			m_A = A;
			m_B = B;
			m_sr = sr;
			m_si = si;
			m_order = order;
			m_modes = modes;
			m_vectors = vectors;
			m_rows_map = rows_map;
			m_cols_map = cols_map;
		}

		//operation
		void SparseNonGen::operation(double* y, double* x) const
		{
			if(m_ido == 1 || m_ido == -1)
			{
				math::Sparse(m_A, m_rows_map, m_cols_map, m_order, m_order).product(y, x);
				memcpy(x, y, m_order * sizeof(double));
				m_S->substitute(y, x ,1);
			}
			if(m_ido == 2)
			{
				math::Sparse(m_B, m_rows_map, m_cols_map, m_order, m_order).product(y, x);
			}
		}
	}
}