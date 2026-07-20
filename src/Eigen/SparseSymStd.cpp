//std
#include <cstring>

//Math
#include "Math/inc/Linear/Sparse.hpp"
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
		SparseSymStd::SparseSymStd(void) : SparseSym{0U, 0U, 0U, nullptr, nullptr}, m_s{nullptr}, m_U{nullptr}, m_A{nullptr}
		{
			return;
		}
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
				dsaupd_(&m_ido, "I", &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
				//compute
				if(m_ido == 99) break;
				if(m_ido == -1 || m_ido == 1) operation(workd + ipntr[1] - 1, workd + ipntr[0] - 1);
			}
			if(info != 0) return false;
			//eigenvalues
			const uint32_t rvec = bool(m_U);
			uint32_t* select = new uint32_t[ncv];
			dseupd_(&rvec, "A", select, m_s, m_U, &n, &m_shift, "I", &n, which, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
			//delete
			delete[] v;
			delete[] resid;
			delete[] workd;
			delete[] workl;
			delete[] select;
			//return
			return info == 0;
		}
		void SparseSymStd::setup(uint32_t order, uint32_t modes, uint32_t vectors, const int32_t* rows_map, const int32_t* cols_map, const double* A, double* s, double* U)
		{
			m_s = s;
			m_U = U;
			m_A = A;
			m_order = order;
			m_modes = modes;
			m_vectors = vectors;
			m_rows_map = rows_map;
			m_cols_map = cols_map;
		}

		//operation
		void SparseSymStd::operation(double* y, double* x) const
		{
			math::Sparse(m_A, m_rows_map, m_cols_map, m_order, m_order).product(y, x);
		}
	}
}