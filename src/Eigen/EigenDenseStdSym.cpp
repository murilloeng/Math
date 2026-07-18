//std
#include <cfloat>
#include <cstring>

//Math
#include "Math/inc/Eigen/EigenDenseStdSym.hpp"

extern "C"
{
	void dsyev_(const char*, const char*, const uint32_t*, double*, const uint32_t*, double*, double*, const int32_t*, int32_t*);
	void dsyevx_(const char*, const char*, const char*, const uint32_t*, double*, const uint32_t*, const double*, const double*, const uint32_t*, const uint32_t*, const double*, uint32_t*, double*, double*, const uint32_t*, double*, const int32_t*, int32_t*, int32_t*, int32_t*);
}

namespace math
{
	namespace eigen
	{
		//constructor
		EigenDenseStdSym::EigenDenseStdSym(uint32_t order, double* A, double* s, double* U) : EigenDenseStd{order, A}, 
			m_s{s}, m_U{U}, m_range{'A'}, m_value_min{-DBL_MAX}, m_value_max{+DBL_MAX}, m_index_min{0}, m_index_max{order}
		{
			return;
		}

		//destructor
		EigenDenseStdSym::~EigenDenseStdSym(void)
		{
			return;
		}

		//data
		char EigenDenseStdSym::range(void) const
		{
			return m_range;
		}
		char EigenDenseStdSym::range(char range)
		{
			return m_range = range;
		}

		uint32_t EigenDenseStdSym::modes(void) const
		{
			return m_modes;
		}

		double EigenDenseStdSym::value_min(void) const
		{
			return m_value_min;
		}
		double EigenDenseStdSym::value_min(double value_min)
		{
			return m_value_min = value_min;
		}

		double EigenDenseStdSym::value_max(void) const
		{
			return m_value_max;
		}
		double EigenDenseStdSym::value_max(double value_max)
		{
			return m_value_max = value_max;
		}

		uint32_t EigenDenseStdSym::index_min(void) const
		{
			return m_index_min;
		}
		uint32_t EigenDenseStdSym::index_min(uint32_t index_min)
		{
			return m_index_min = index_min;
		}

		uint32_t EigenDenseStdSym::index_max(void) const
		{
			return m_index_max;
		}
		uint32_t EigenDenseStdSym::index_max(uint32_t index_max)
		{
			return m_index_max = index_max;
		}

		//compute
		bool EigenDenseStdSym::compute(void)
		{
			return m_range == 'A' ? compute_full() : compute_partial();
		}
		bool EigenDenseStdSym::compute_full(void)
		{
			//data
			double query;
			int32_t info;
			int32_t lwork = -1;
			const char uplo = 'U';
			const char jobz = m_U ? 'V' : 'N';
			//query
			dsyev_(&jobz, &uplo, &m_order, m_A, &m_order, m_s, &query, &lwork, &info);
			//compute
			lwork = int32_t(query);
			double* work = new double[lwork];
			dsyev_(&jobz, &uplo, &m_order, m_A, &m_order, m_s, work, &lwork, &info);
			//eigenvectors
			if(m_U) memcpy(m_U, m_A, m_order * m_order * sizeof(double));
			//delete
			delete[] work;
			//return
			return info == 0;
		}
		bool EigenDenseStdSym::compute_partial(void)
		{
			//data
			double query;
			int32_t info;
			int32_t lwork = -1;
			const char uplo = 'U';
			const char jobz = m_U ? 'V' : 'N';
			//query
			const double abstol = 0;
			const uint32_t* n = &m_order;
			const double v1 = m_value_min;
			const double v2 = m_value_max;
			const uint32_t i1 = m_index_min + 1;
			const uint32_t i2 = m_index_max + 1;
			int32_t* ifail = new int32_t[m_order];
			int32_t* iwork = new int32_t[5 * m_order];
			dsyevx_(&jobz, &m_range, &uplo, n, m_A, n, &v1, &v2, &i1, &i2, &abstol, &m_modes, m_s, m_U, n, &query, &lwork, iwork, ifail, &info);
			//compute
			lwork = int32_t(query);
			double* work = new double[lwork];
			dsyevx_(&jobz, &m_range, &uplo, n, m_A, n, &v1, &v2, &i1, &i2, &abstol, &m_modes, m_s, m_U, n, work, &lwork, iwork, ifail, &info);
			//delete
			delete[] work;
			delete[] ifail;
			delete[] iwork;
			//return
			return info == 0;
		}
	}
}