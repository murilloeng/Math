//std
#include <cfloat>
#include <cstring>

//Math
#include "Math/inc/Eigen/DenseNonGen.hpp"

extern "C"
{
	void dggev_(const char*, const char*, const uint32_t*, double*, const uint32_t*, double*, const uint32_t*, double*, double*, double*, double*, const uint32_t*, double*, const uint32_t*, double*, const int32_t*, int32_t*);
}

namespace math
{
	namespace eigen
	{
		//constructor
		DenseNonGen::DenseNonGen(uint32_t order, double* A, double* B, double* sr, double* si, double* U, double* V) : 
			DenseNon{order}, m_A{A}, m_B{B}, m_U{U}, m_V{V}, m_sr{sr}, m_si{si}
		{
			return;
		}

		//destructor
		DenseNonGen::~DenseNonGen(void)
		{
			return;
		}

		//compute
		bool DenseNonGen::compute(void)
		{
			//data
			double query;
			int32_t info;
			int32_t lwork = -1;
			const char jobvl = m_V ? 'V' : 'N';
			const char jobvr = m_U ? 'V' : 'N';
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
	}
}