//std
#include <cmath>
#include <cstring>

//Math
#include "Math/inc/Linear/SVD.hpp"
#include "Math/inc/Miscellaneous/util.hpp"

extern "C"
{
	void dgesvd_(const char*, const char*, const uint32_t*, const uint32_t*, double*, const uint32_t*, double*, double*, const uint32_t*, double*, const uint32_t*, double*, const int32_t*, int32_t*);
}

namespace math
{
	//constructor
	SVD::SVD(const double* A, uint32_t rows, uint32_t cols) : 
		m_own{true}, m_modes{true}, m_rows{rows}, m_cols{cols}, m_A{A}, m_s{nullptr}, m_U{nullptr}, m_V{nullptr}
	{
		return;
	}
	SVD::SVD(const double* A, uint32_t rows, uint32_t cols, double* s, double* U, double* V) : 
		m_own{false}, m_modes{true}, m_rows{rows}, m_cols{cols}, m_A{A}, m_s{s}, m_U{U}, m_V{V}
	{
		return;
	}

	//destructor
	SVD::~SVD(void)
	{
		if(m_own)
		{
			delete[] m_s;
			delete[] m_U;
			delete[] m_V;
		}
	}

	//data
	bool SVD::modes(void) const
	{
		return m_modes;
	}
	bool SVD::modes(bool modes)
	{
		return m_modes = modes;
	}

	uint32_t SVD::rows(void) const
	{
		return m_rows;
	}
	uint32_t SVD::cols(void) const
	{
		return m_cols;
	}
	const double* SVD::data(void) const
	{
		return m_A;
	}
	const double* SVD::singular_values(void) const
	{
		return m_s;
	}
	const double* SVD::singular_vectors(uint32_t index) const
	{
		return index == 0 ? m_U : m_V;
	}

	//compute
	bool SVD::compute(void)
	{
		//data
		double query;
		int32_t info;
		int32_t lwork = -1;
		const char jobu = m_modes ? 'A' : 'N';
		const char jobv = m_modes ? 'A' : 'N';
		double* A = new double[m_rows * m_cols];
		memcpy(A, m_A, m_rows * m_cols * sizeof(double));
		//query
		cleanup();
		allocate();
		dgesvd_(&jobu, &jobv, &m_rows, &m_cols, A, &m_rows, m_s, m_U, &m_rows, m_V, &m_cols, &query, &lwork, &info);
		//compute
		lwork = int32_t(query);
		double* work = new double[lwork];
		dgesvd_(&jobu, &jobv, &m_rows, &m_cols, A, &m_rows, m_s, m_U, &m_rows, m_V, &m_cols, work, &lwork, &info);
		//transpose
		for(uint32_t i = 0; i < m_cols; i++)
		{
			for(uint32_t j = i + 1; j < m_cols; j++)
			{
				swap(m_V[i + m_cols * j], m_V[j + m_cols * i]);
			}
		}
		//delete
		delete[] A;
		delete[] work;
		//return
		return info == 0;
	}

	//setup
	void SVD::cleanup(void)
	{
		if(m_own)
		{
			delete[] m_s;
			delete[] m_U;
			delete[] m_V;
		}
	}
	void SVD::allocate(void)
	{
		if(m_own)
		{
			m_U = new double[m_rows * m_rows];
			m_V = new double[m_cols * m_cols];
			m_s = new double[std::min(m_rows, m_cols)];
		}
	}
}