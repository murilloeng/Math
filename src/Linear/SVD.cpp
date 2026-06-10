//std
#include <cmath>
#include <cstring>

//Math
#include "Math/inc/Linear/SVD.hpp"

extern "C"
{
	void dgesvd_(const char*, const char*, const uint32_t*, const uint32_t*, double*, const uint32_t*, double*, double*, const uint32_t*, double*, const uint32_t*, double*, const int32_t*, int32_t*);
}

namespace math
{
	//constructor
	SVD::SVD(void) : 
		m_modes{true}, m_rows{0}, m_cols{0}, m_data{nullptr}, 
		m_singular_values{nullptr}, m_singular_modes{nullptr, nullptr}
	{
		return;
	}
	
	//destructor
	SVD::~SVD(void)
	{
		delete[] m_singular_values;
		delete[] m_singular_modes[0];
		delete[] m_singular_modes[1];
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
	uint32_t SVD::rows(uint32_t rows)
	{
		return m_rows = rows;
	}

	uint32_t SVD::cols(void) const
	{
		return m_cols;
	}
	uint32_t SVD::cols(uint32_t cols)
	{
		return m_cols = cols;
	}

	const double* SVD::data(void) const
	{
		return m_data;
	}
	const double* SVD::data(const double* data)
	{
		return m_data = data;
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
		memcpy(A, m_data, m_rows * m_cols * sizeof(double));
		//query
		double* s = m_singular_values;
		double* U = m_singular_modes[0];
		double* V = m_singular_modes[1];
		dgesvd_(&jobu, &jobv, &m_rows, &m_cols, A, &m_rows, s, U, &m_rows, V, &m_cols, &query, &lwork, &info);
		//compute
		lwork = int32_t(query);
		double* work = new double[lwork];
		dgesvd_(&jobu, &jobv, &m_rows, &m_cols, A, &m_rows, s, U, &m_rows, V, &m_cols, work, &lwork, &info);
		//delete
		delete[] A;
		delete[] work;
		//return
		return info == 0;
	}

	//singular values
	const double* SVD::singular_values(void) const
	{
		return m_singular_values;
	}
	double SVD::singular_value(uint32_t index) const
	{
		return m_singular_values[index];
	}

	//singular modes
	const double* SVD::singular_modes(uint32_t type) const
	{
		return m_singular_modes[type];
	}
	const double* SVD::singular_modes(uint32_t type, uint32_t index) const
	{
		return m_singular_modes[type] + index * m_rows;
	}

	//setup
	void SVD::cleanup(void)
	{
		delete[] m_singular_values;
		delete[] m_singular_modes[0];
		delete[] m_singular_modes[1];
	}
	void SVD::allocate(void)
	{
		//data
		const uint32_t n = std::min(m_rows, m_cols);
		//allocate
		m_singular_values = new double[n];
		m_singular_modes[0] = new double[m_rows * m_rows];
		m_singular_modes[1] = new double[m_cols * m_cols];
	}
}