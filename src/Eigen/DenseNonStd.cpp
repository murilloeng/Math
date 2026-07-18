//std
#include <cfloat>
#include <cstring>

//Math
#include "Math/inc/Eigen/DenseNonStd.hpp"

extern "C"
{
	void dgeev_(const char*, const char*, const uint32_t*, double*, const uint32_t*, double*, double*, double*, const uint32_t*, double*, const uint32_t*, double*, int32_t*, int32_t*);
}

namespace math
{
	namespace eigen
	{
		//constructor
		DenseNonStd::DenseNonStd(uint32_t order, double* A, double* sr, double* si, double* U, double* V) : 
			DenseNon{order}, m_A{A}, m_U{U}, m_V{V}, m_sr{sr}, m_si{si}
		{
			return;
		}

		//destructor
		DenseNonStd::~DenseNonStd(void)
		{
			return;
		}

		//compute
		bool DenseNonStd::compute(void)
		{
			//data
			double query;
			int32_t info;
			int32_t lwork = -1;
			const char jobvl = m_V ? 'V' : 'N';
			const char jobvr = m_U ? 'V' : 'N';
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
	}
}