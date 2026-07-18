//std
#include <cstring>

//Math
#include "Math/inc/Eigen/SparseSymStd.hpp"

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
		SparseSymStd::SparseSymStd(uint32_t order, uint32_t modes, uint32_t vectors, const int32_t* rows_map, const int32_t* cols_map, const double* A, double* s, double* U) : SparseSym{order, modes, vectors, rows_map, cols_map}, m_s{s}, m_U{U}, m_A{A}
		{
			return;
		}

		//destructor
		SparseSymStd::~SparseSymStd(void)
		{
			return;
		}

		//compute
		bool SparseSymStd::compute(void)
		{
			//data
			const uint32_t n = m_order;
			const uint32_t nev = m_modes;
			const uint32_t ncv = m_vectors;
			const uint32_t lworkl = (ncv + 8) * ncv;
			//data
			const double tol = 1.00e-10;
			double* resid = new double[n];
			double* v = new double[n * ncv];
			double* workd = new double[3 * n];
			double* workl = new double[lworkl];
			int32_t ido = 0, info = 0, iparam[11], ipntr[11];
			//setup
			iparam[0] = 1;
			iparam[6] = 1;
			iparam[2] = 300;
			//Arnoldi
			while(true)
			{
				//decomposition
				dsaupd_(&ido, "I", &n, "SA", &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
				//compute
				if(ido == 99) break;
				if(ido == -1 || ido == 1) sparse_product(workd + ipntr[1] - 1, workd + ipntr[0] - 1);
			}
			if(info != 0) return false;
			//eigenvalues
			const uint32_t rvec = true;
			uint32_t* select = new uint32_t[ncv];
			dseupd_(&rvec, "A", select, m_s, m_U, &n, nullptr, "I", &n, "SA", &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
			//delete
			delete[] v;
			delete[] resid;
			delete[] workd;
			delete[] workl;
			delete[] select;
			//return
			return info == 0;
		}

		//operation
		void SparseSymStd::sparse_product(double* y, const double* x) const
		{
			memset(y, 0, m_order * sizeof(double));
			for(uint32_t i = 0; i < m_order; i++)
			{
				for(int32_t j = m_cols_map[i]; j < m_cols_map[i + 1]; j++)
				{
					y[m_rows_map[j]] += m_A[j] * x[i];
				}
			}
		}
	}
}