//std
#include <cfloat>
#include <cstdio>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Linear/Eigen.hpp"

extern "C"
{
	void dsyev_(const char*, const char*, const uint32_t*, double*, const uint32_t*, double*, double*, const int32_t*, int32_t*);
	void dsygv_(const uint32_t*, const char*, const char*, const uint32_t*, double*, const uint32_t*, double*, const uint32_t*, double*, double*, const int32_t*, int32_t*);
	void dgeev_(const char*, const char*, const uint32_t*, double*, const uint32_t*, double*, double*, double*, const uint32_t*, double*, const uint32_t*, double*, int32_t*, int32_t*);
	void dggev_(const char*, const char*, const uint32_t*, double*, const uint32_t*, double*, const uint32_t*, double*, double*, double*, double*, const uint32_t*, double*, const uint32_t*, double*, const int32_t*, int32_t*);
	void dsyevx_(const char*, const char*, const char*, const uint32_t*, double*, const uint32_t*, const double*, const double*, const uint32_t*, const uint32_t*, const double*, uint32_t*, double*, double*, const uint32_t*, double*, const int32_t*, int32_t*, int32_t*, int32_t*);
	void dsygvx_(const uint32_t*, const char*, const char*, const char*, const uint32_t*, double*, const uint32_t*, double*, const uint32_t*, const double*, const double*, const uint32_t*, const uint32_t*, const double*, uint32_t*, double*, double*, const uint32_t*, double*, const int32_t*, int32_t*, int32_t*, int32_t*);
	void dsaupd_(int32_t*, const char*, const uint32_t*, const char*, const uint32_t*, const double*, double*, const uint32_t*, double*, const uint32_t*, int32_t*, int32_t*, double*, double*, const uint32_t*, int32_t*);
	void dseupd_(const uint32_t*, const char*, uint32_t*, double*, double*, const uint32_t*, const double*, const char*, const uint32_t*, const char*, const uint32_t*, const double*, double*, const uint32_t*, double*, const uint32_t*, int32_t*, int32_t*, double*, double*, const uint32_t*, int32_t*);
}

namespace math
{
	//constructor
	Eigen::Eigen(double* A, uint32_t order, double* s, double* U) :
		m_full{true}, m_dense{true}, m_symmetric{true}, m_order{order}, m_index_min{0}, m_index_max{order - 1},
		m_rows_map{nullptr}, m_cols_map{nullptr}, m_A{A}, m_B{nullptr}, m_sr{s}, m_si{nullptr}, m_U{U}, m_V{nullptr}
	{
		return;
	}
	Eigen::Eigen(double* A, uint32_t order, double* sr, double* si, double* U, double* V) :
		m_full{true}, m_dense{true}, m_symmetric{false}, m_order{order}, m_index_min{0}, m_index_max{order - 1},
		m_rows_map{nullptr}, m_cols_map{nullptr}, m_A{A}, m_B{nullptr}, m_sr{sr}, m_si{si}, m_U{U}, m_V{V}
	{
		return;
	}
	Eigen::Eigen(double* A, uint32_t order, double* s, double* U, uint32_t index_min, uint32_t index_max) :
		m_full{false}, m_dense{true}, m_symmetric{true}, m_order{order}, m_index_min{index_min}, m_index_max{index_max},
		m_rows_map{nullptr}, m_cols_map{nullptr}, m_A{A}, m_B{nullptr}, m_sr{s}, m_si{nullptr}, m_U{U}, m_V{nullptr}
	{
		return;
	}

	Eigen::Eigen(double* A, double* B, uint32_t order, double* s, double* U) :
		m_full{true}, m_dense{true}, m_symmetric{true}, m_order{order}, m_index_min{0}, m_index_max{order - 1},
		m_rows_map{nullptr}, m_cols_map{nullptr}, m_A{A}, m_B{B}, m_sr{s}, m_si{nullptr}, m_U{U}, m_V{nullptr}
	{
		return;
	}
	Eigen::Eigen(double* A, double* B, uint32_t order, double* sr, double* si, double* U, double* V) :
		m_full{true}, m_dense{true}, m_symmetric{false}, m_order{order}, m_index_min{0}, m_index_max{order - 1},
		m_rows_map{nullptr}, m_cols_map{nullptr}, m_A{A}, m_B{B}, m_sr{sr}, m_si{si}, m_U{U}, m_V{V}
	{
		return;
	}
	Eigen::Eigen(double* A, double* B, uint32_t order, double* s, double* U, uint32_t index_min, uint32_t index_max) :
		m_full{true}, m_dense{true}, m_symmetric{true}, m_order{order},m_index_min{index_min}, m_index_max{index_max},
		m_rows_map{nullptr}, m_cols_map{nullptr}, m_A{A}, m_B{B}, m_sr{s}, m_si{nullptr}, m_U{U}, m_V{nullptr}
	{
		return;
	}

	Eigen::Eigen(double* A, uint32_t order, const int32_t* rows_map, const int32_t* cols_map, double* s, double* U, uint32_t nev, uint32_t ncv) : 
		m_full{false}, m_dense{false}, m_symmetric{true}, m_order{order}, m_index_min{0}, m_index_max{0}, m_nev{nev}, m_ncv{ncv},
		m_rows_map{rows_map}, m_cols_map{cols_map}, m_A{A}, m_B{nullptr}, m_sr{s}, m_si{nullptr}, m_U{U}, m_V{nullptr}
	{
		return;
	}
	Eigen::Eigen(double* A, uint32_t order, const int32_t* rows_map, const int32_t* cols_map, double* sr, double* si, double* U, double* V, uint32_t nev, uint32_t ncv) : 
		m_full{false}, m_dense{false}, m_symmetric{false}, m_order{order}, m_index_min{0}, m_index_max{0}, m_nev{nev}, m_ncv{ncv},
		m_rows_map{rows_map}, m_cols_map{cols_map}, m_A{A}, m_B{nullptr}, m_sr{sr}, m_si{si}, m_U{U}, m_V{V}
	{
		return;
	}

	Eigen::Eigen(double* A, double* B, uint32_t order, const int32_t* rows_map, const int32_t* cols_map, double* s, double* U, uint32_t nev, uint32_t ncv) : 
		m_full{false}, m_dense{false}, m_symmetric{true}, m_order{order}, m_index_min{0}, m_index_max{0}, m_nev{nev}, m_ncv{ncv},
		m_rows_map{rows_map}, m_cols_map{cols_map}, m_A{A}, m_B{B}, m_sr{s}, m_si{nullptr}, m_U{U}, m_V{nullptr}
	{
		return;
	}
	Eigen::Eigen(double* A, double* B, uint32_t order, const int32_t* rows_map, const int32_t* cols_map, double* sr, double* si, double* U, double* V, uint32_t nev, uint32_t ncv) : 
		m_full{false}, m_dense{false}, m_symmetric{false}, m_order{order}, m_index_min{0}, m_index_max{0}, m_nev{nev}, m_ncv{ncv},
		m_rows_map{rows_map}, m_cols_map{cols_map}, m_A{A}, m_B{B}, m_sr{sr}, m_si{si}, m_U{U}, m_V{V}
	{
		return;
	}

	//destructor
	Eigen::~Eigen(void)
	{
		return;
	}

	//data
	bool Eigen::full(void) const
	{
		return m_full;
	}

	bool Eigen::modes(void) const
	{
		return m_symmetric ? bool(m_U) : m_U && m_V;
	}

	bool Eigen::symmetric(void) const
	{
		return m_symmetric;
	}

	uint32_t Eigen::order(void) const
	{
		return m_order;
	}

	double* Eigen::data(uint32_t index) const
	{
		return index == 0 ? m_A : m_B;
	}

	uint32_t Eigen::index_min(void) const
	{
		return m_index_min;
	}
	uint32_t Eigen::index_max(void) const
	{
		return m_index_max;
	}

	//compute
	bool Eigen::compute(void)
	{
		//data
		const bool std = m_B == nullptr;
		bool(Eigen::*methods[])(void) = {
			&Eigen::compute_dense_symmetric_std_full,
			&Eigen::compute_dense_symmetric_gen_full,
			&Eigen::compute_dense_symmetric_std_partial,
			&Eigen::compute_dense_symmetric_gen_partial,
			&Eigen::compute_dense_non_symmetric_std_full,
			&Eigen::compute_dense_non_symmetric_gen_full,
			&Eigen::compute_sparse_symmetric_std_partial,
			&Eigen::compute_sparse_symmetric_gen_partial,
			&Eigen::compute_sparse_non_symmetric_std_partial,
			&Eigen::compute_sparse_non_symmetric_gen_partial
		};
		const uint32_t index = 
			 m_dense &&  m_symmetric &&  std &&  m_full ? 0 :
			 m_dense &&  m_symmetric && !std &&  m_full ? 1 :
			 m_dense &&  m_symmetric &&  std && !m_full ? 2 :
			 m_dense &&  m_symmetric && !std && !m_full ? 3 :
			 m_dense && !m_symmetric &&  std &&  m_full ? 4 :
			 m_dense && !m_symmetric && !std &&  m_full ? 5 :
			!m_dense &&  m_symmetric &&  std && !m_full ? 6 :
			!m_dense &&  m_symmetric && !std && !m_full ? 7 :
			!m_dense && !m_symmetric &&  std && !m_full ? 8 :
			!m_dense && !m_symmetric && !std && !m_full ? 9 : 0;
		//compute
		return (this->*methods[index])();
	}

	//eigenvalues
	const double* Eigen::eigenvalues(uint32_t type) const
	{
		return type == 0 ? m_sr : m_si;
	}
	double Eigen::eigenvalue(uint32_t type, uint32_t index) const
	{
		return (type == 0 ? m_sr : m_si)[index];
	}

	//eigenvectors
	const double* Eigen::eigenvectors(uint32_t type) const
	{
		return type == 0 ? m_U : m_V;
	}
	const double* Eigen::eigenvector(uint32_t type, uint32_t index) const
	{
		return (type == 0 ? m_U : m_V) + index * m_order;
	}

	//compute
	bool Eigen::compute_dense_symmetric_std_full(void)
	{
		//data
		double query;
		int32_t info;
		int32_t lwork = -1;
		const char uplo = 'U';
		const char jobz = m_U ? 'V' : 'N';
		memcpy(m_U, m_A, m_order * m_order * sizeof(double));
		//query
		dsyev_(&jobz, &uplo, &m_order, m_U, &m_order, m_sr, &query, &lwork, &info);
		//compute
		lwork = int32_t(query);
		double* work = new double[lwork];
		dsyev_(&jobz, &uplo, &m_order, m_U, &m_order, m_sr, work, &lwork, &info);
		//delete
		delete[] work;
		//return
		return info == 0;
	}
	bool Eigen::compute_dense_symmetric_gen_full(void)
	{
		//data
		double query;
		int32_t info;
		int32_t lwork = -1;
		uint32_t itype = 1;
		const char uplo = 'U';
		const char jobz = m_U ? 'V' : 'N';
		memcpy(m_U, m_A, m_order * m_order * sizeof(double));
		//query
		dsygv_(&itype, &jobz, &uplo, &m_order, m_U, &m_order, m_B, &m_order, m_sr, &query, &lwork, &info);
		//compute
		lwork = int32_t(query);
		double* work = new double[lwork];
		dsygv_(&itype, &jobz, &uplo, &m_order, m_U, &m_order, m_B, &m_order, m_sr, work, &lwork, &info);
		//delete
		delete[] work;
		//return
		return info == 0;
	}
	bool Eigen::compute_dense_symmetric_std_partial(void)
	{
		//data
		uint32_t m;
		double query;
		int32_t info;
		int32_t lwork = -1;
		const char uplo = 'U';
		const char range = 'I';
		const char jobz = m_U ? 'V' : 'N';
		//query
		const double abstol = 0;
		const double* v1 = nullptr;
		const double* v2 = nullptr;
		const uint32_t* n = &m_order;
		const uint32_t i1 = m_index_min + 1;
		const uint32_t i2 = m_index_max + 1;
		int32_t* ifail = new int32_t[m_order];
		int32_t* iwork = new int32_t[5 * m_order];
		dsyevx_(&jobz, &range, &uplo, n, m_A, n, v1, v2, &i1, &i2, &abstol, &m, m_sr, m_U, n, &query, &lwork, iwork, ifail, &info);
		//compute
		lwork = int32_t(query);
		double* work = new double[lwork];
		dsyevx_(&jobz, &range, &uplo, n, m_A, n, v1, v2, &i1, &i2, &abstol, &m, m_sr, m_U, n, work, &lwork, iwork, ifail, &info);
		//delete
		delete[] work;
		delete[] ifail;
		delete[] iwork;
		//return
		return info == 0;
	}
	bool Eigen::compute_dense_symmetric_gen_partial(void)
	{
		//data
		uint32_t m;
		double query;
		int32_t info;
		int32_t lwork = -1;
		const char uplo = 'U';
		const char range = 'I';
		const char jobz = m_U ? 'V' : 'N';
		//query
		const double abstol = 0;
		const uint32_t itype = 1;
		const double* v1 = nullptr;
		const double* v2 = nullptr;
		const uint32_t* n = &m_order;
		const uint32_t i1 = m_index_min + 1;
		const uint32_t i2 = m_index_max + 1;
		int32_t* ifail = new int32_t[m_order];
		int32_t* iwork = new int32_t[5 * m_order];
		dsygvx_(&itype, &jobz, &range, &uplo, n, m_A, n, m_B, n, v1, v2, &i1, &i2, &abstol, &m, m_sr, m_U, n, &query, &lwork, iwork, ifail, &info);
		//compute
		lwork = int32_t(query);
		double* work = new double[lwork];
		dsygvx_(&itype, &jobz, &range, &uplo, n, m_A, n, m_B, n, v1, v2, &i1, &i2, &abstol, &m, m_sr, m_U, n, work, &lwork, iwork, ifail, &info);
		//delete
		delete[] work;
		delete[] ifail;
		delete[] iwork;
		//return
		return info == 0;
	}
	bool Eigen::compute_dense_non_symmetric_std_full(void)
	{
		//data
		double query;
		int32_t info;
		int32_t lwork = -1;
		const char jobvl = m_U && m_V ? 'V' : 'N';
		const char jobvr = m_U && m_V ? 'V' : 'N';
		//query
		const uint32_t* n = &m_order;
		dgeev_(&jobvl, &jobvr, n, m_A, n, m_sr, m_si, m_V, n, m_U, n, &query, &lwork, &info);
		//compute
		lwork = int32_t(query);
		double* work = new double[lwork];
		dgeev_(&jobvl, &jobvr, n, m_A, n, m_sr, m_si, m_V, n, m_U, n, work, &lwork, &info);
		//delete
		delete[] work;
		//return
		return info == 0;
	}
	bool Eigen::compute_dense_non_symmetric_gen_full(void)
	{
		//data
		double query;
		int32_t info;
		int32_t lwork = -1;
		const char jobvl = m_U && m_V ? 'V' : 'N';
		const char jobvr = m_U && m_V ? 'V' : 'N';
		//query
		const uint32_t* n = &m_order;
		double* b = new double[m_order];
		dggev_(&jobvl, &jobvr, n, m_A, n, m_B, n, m_sr, m_si, b, m_V, n, m_U, n, &query, &lwork, &info);
		//compute
		lwork = int32_t(query);
		double* work = new double[lwork];
		dggev_(&jobvl, &jobvr, n, m_A, n, m_B, n, m_sr, m_si, b, m_V, n, m_U, n, work, &lwork, &info);
		//division
		for(uint32_t i = 0; i < m_order; i++)
		{
			m_sr[i] /= b[i];
			m_si[i] /= b[i];
		}
		//delete
		delete[] b;
		delete[] work;
		//return
		return info == 0;
	}

	bool Eigen::compute_sparse_symmetric_std_partial(void)
	{
		//data
		const uint32_t n = m_order;
		const uint32_t lworkl = (m_ncv + 8) * m_ncv;
		//data
		const double tol = 1.00e-10;
		double* resid = new double[n];
		double* v = new double[n * m_ncv];
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
			dsaupd_(&ido, "I", &n, "SA", &m_nev, &tol, resid, &m_ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
			//compute
			if(ido == 99) break;
			if(ido == -1 || ido == 1) sparse_product(workd + ipntr[1] - 1, workd + ipntr[0] - 1);
		}
		if(info != 0) return false;
		//eigenvalues
		const uint32_t rvec = true;
		uint32_t* select = new uint32_t[m_ncv];
		dseupd_(&rvec, "A", select, m_sr, m_U, &n, nullptr, "I", &n, "SA", &m_nev, &tol, resid, &m_ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, &info);
		//delete
		delete[] v;
		delete[] resid;
		delete[] workd;
		delete[] workl;
		delete[] select;
		//return
		return info == 0;
	}
	bool Eigen::compute_sparse_symmetric_gen_partial(void)
	{
		return false;
	}
	bool Eigen::compute_sparse_non_symmetric_std_partial(void)
	{
		return false;
	}
	bool Eigen::compute_sparse_non_symmetric_gen_partial(void)
	{
		return false;
	}

	//sparse
	void Eigen::sparse_product(double* y, const double* x) const
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