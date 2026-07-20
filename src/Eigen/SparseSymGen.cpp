//std
#include <cstdio>
#include <cstring>

//Math
#include "Math/inc/Linear/Sparse.hpp"
#include "Math/inc/Eigen/SparseSymGen.hpp"

extern "C"
{
	void dsaupd_(int32_t*, const char*, const uint32_t*, const char*, const uint32_t*, const double*, double*, const uint32_t*, double*, const uint32_t*, int32_t*, int32_t*, double*, double*, const uint32_t*, int32_t*);
	void dseupd_(const uint32_t*, const char*, uint32_t*, double*, double*, const uint32_t*, const double*, const char*, const uint32_t*, const char*, const uint32_t*, const double*, double*, const uint32_t*, double*, const uint32_t*, int32_t*, int32_t*, double*, double*, const uint32_t*, int32_t*);
}

namespace math
{
	namespace eigen
	{
		//constructor
		SparseSymGen::SparseSymGen(uint32_t order, uint32_t modes, uint32_t vectors, const int32_t* rows_map, const int32_t* cols_map, const double* A, const double* B, double* s, double* U) : SparseSym{order, modes, vectors, rows_map, cols_map}, m_s{s}, m_U{U}, m_A{A}, m_B{B}
		{
			m_mode = 2;
		}

		//destructor
		SparseSymGen::~SparseSymGen(void)
		{
			return;
		}

		//compute
		bool SparseSymGen::compute(void)
		{
			//data
			const uint32_t n = m_order;
			const uint32_t nev = m_modes;
			const uint32_t ncv = m_vectors;
			const uint32_t lworkl = (ncv + 8) * ncv;
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
				dsaupd_(&m_ido, "G", &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
				//compute
				if(m_ido == 99) break;
				if(m_ido == -1 || m_ido == 1 || m_ido == 2) operation(workd + ipntr[1] - 1, workd + ipntr[0] - 1);
			}
			if(info != 0) return false;
			//eigenvalues
			const uint32_t rvec = bool(m_U);
			uint32_t* select = new uint32_t[ncv];
			dseupd_(&rvec, "A", select, m_s, m_U, &n, nullptr, "G", &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
			//delete
			delete[] v;
			delete[] resid;
			delete[] workd;
			delete[] workl;
			delete[] select;
			m_S->release();
			//return
			return info == 0;
		}

		//operation
		void SparseSymGen::operation(double* y, double* x) const
		{
			if(m_ido == 1 || m_ido == -1)
			{
				math::Sparse(m_A, m_rows_map, m_cols_map, m_order, m_order).product(y, x);
				memcpy(x, y, m_order * sizeof(double));
				// !math::Sparse(m_B, m_rows_map, m_cols_map, m_order, m_order).solve(y, x, 1);
				m_S->substitute(y, x ,1);
			}
			if(m_ido == 2)
			{
				math::Sparse(m_B, m_rows_map, m_cols_map, m_order, m_order).product(y, x);
			}
		}
	}
}