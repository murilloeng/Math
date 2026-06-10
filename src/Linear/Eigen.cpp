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
}

namespace math
{
	//constructor
	Eigen::Eigen(void) : 
		m_type{Type::Full}, m_symmetric{true}, m_order{0},
		m_value_min{-DBL_MAX}, m_value_max{+DBL_MAX}, m_index_min{0}, m_index_max{UINT32_MAX},
		m_data{nullptr, nullptr}, m_eigenvalues{nullptr, nullptr}, m_eigenvectors{nullptr, nullptr}, m_eigenvectors_computation{true, true}
	{
		return;
	}
	
	//destructor
	Eigen::~Eigen(void)
	{
		return;
	}

	//data
	Eigen::Type Eigen::type(void) const
	{
		return m_type;
	}
	Eigen::Type Eigen::type(Type type)
	{
		return m_type = type;
	}

	bool Eigen::symmetric(void) const
	{
		return m_symmetric;
	}
	bool Eigen::symmetric(bool symmetry)
	{
		return m_symmetric = symmetry;
	}

	uint32_t Eigen::order(void) const
	{
		return m_order;
	}
	uint32_t Eigen::order(uint32_t order)
	{
		return m_order = order;
	}

	double Eigen::value_min(void) const
	{
		return m_value_min;
	}
	double Eigen::value_min(double value_min)
	{
		return m_value_min = value_min;
	}

	double Eigen::value_max(void) const
	{
		return m_value_max;
	}
	double Eigen::value_max(double value_max)
	{
		return m_value_max = value_max;
	}

	uint32_t Eigen::index_min(void) const
	{
		return m_index_min;
	}
	uint32_t Eigen::index_min(uint32_t index_min)
	{
		return m_index_min = index_min;
	}

	uint32_t Eigen::index_max(void) const
	{
		return m_index_max;
	}
	uint32_t Eigen::index_max(uint32_t index_max)
	{
		return m_index_max = index_max;
	}

	const double* Eigen::data(uint32_t index) const
	{
		return m_data[index];
	}
	const double* Eigen::data(uint32_t index, const double* data)
	{
		return m_data[index] = data;
	}

	bool Eigen::eigenvectors_computation(uint32_t index) const
	{
		return m_eigenvectors_computation[index];
	}
	bool Eigen::eigenvectors_computation(uint32_t index, bool eigenvectors_computation)
	{
		return m_eigenvectors_computation[index] = eigenvectors_computation;
	}

	//compute
	bool Eigen::compute(void)
	{
		//data
		const bool type = m_type != Type::Full;
		const bool data = m_data[1] != nullptr;
		bool(Eigen::*methods[])(void) = {
			&Eigen::compute_non_symmetric_std, &Eigen::compute_non_symmetric_gen,
			&Eigen::compute_symmetric_std_full, &Eigen::compute_symmetric_std_partial,
			&Eigen::compute_symmetric_gen_full, &Eigen::compute_symmetric_gen_partial
		};
		const uint32_t index = !m_symmetric ? data : 2 + 2 * data + type;
		//check
		if(!m_symmetric && m_type != Type::Full)
		{
			throw std::runtime_error("Error: Selected eigenvalues are not supported for non-symmetric matrices!");
		}
		//compute
		cleanup();
		allocate();
		return (this->*methods[index])();
	}

	//modes
	uint32_t Eigen::modes(void) const
	{
		return m_modes;
	}

	//eigenvalues
	const double* Eigen::eigenvalues(uint32_t part) const
	{
		return m_eigenvalues[part];
	}
	double Eigen::eigenvalue(uint32_t part, uint32_t index) const
	{
		return m_eigenvalues[part][index];
	}

	//eigenvectors
	const double* Eigen::eigenvectors(uint32_t type) const
	{
		return m_eigenvectors[type];
	}
	const double* Eigen::eigenvector(uint32_t type, uint32_t index) const
	{
		return m_eigenvectors[type] + index * m_order;
	}

	//setup
	void Eigen::cleanup(void)
	{
		delete[] m_eigenvalues[0];
		delete[] m_eigenvalues[1];
		delete[] m_eigenvectors[0];
		delete[] m_eigenvectors[1];
	}
	void Eigen::allocate(void)
	{
		m_eigenvalues[0] = new double[m_order];
		m_eigenvectors[0] = new double[m_order * m_order];
		if(!m_symmetric) m_eigenvalues[1] = new double[m_order];
		if(!m_symmetric) m_eigenvectors[1] = new double[m_order * m_order];
	}

	//compute
	bool Eigen::compute_non_symmetric_std(void)
	{
		//data
		double query;
		int32_t info;
		int32_t lwork = -1;
		const char jobvl = !m_eigenvectors_computation[1] ? 'N' : 'V';
		const char jobvr = !m_eigenvectors_computation[0] ? 'N' : 'V';
		//query
		const uint32_t* n = &m_order;
		double* wr = m_eigenvalues[0];
		double* wi = m_eigenvalues[1];
		double* Z1 = m_eigenvectors[0];
		double* Z2 = m_eigenvectors[1];
		double* A = new double[m_order * m_order];
		memcpy(A, m_data[0], m_order * m_order * sizeof(double));
		dgeev_(&jobvl, &jobvr, n, A, n, wr, wi, Z2, n, Z1, n, &query, &lwork, &info);
		//compute
		lwork = int32_t(query);
		double* work = new double[lwork];
		dgeev_(&jobvl, &jobvr, n, A, n, wr, wi, Z2, n, Z1, n, work, &lwork, &info);
		//delete
		delete[] A;
		delete[] work;
		//return
		return info == 0;
	}
	bool Eigen::compute_non_symmetric_gen(void)
	{
		//data
		double query;
		int32_t info;
		int32_t lwork = -1;
		const char jobvl = !m_eigenvectors_computation[1] ? 'N' : 'V';
		const char jobvr = !m_eigenvectors_computation[0] ? 'N' : 'V';
		//query
		const uint32_t* n = &m_order;
		double* wr = m_eigenvalues[0];
		double* wi = m_eigenvalues[1];
		double* Z1 = m_eigenvectors[0];
		double* Z2 = m_eigenvectors[1];
		double* b = new double[m_order];
		double* A = new double[m_order * m_order];
		double* B = new double[m_order * m_order];
		memcpy(A, m_data[0], m_order * m_order * sizeof(double));
		memcpy(B, m_data[1], m_order * m_order * sizeof(double));
		dggev_(&jobvl, &jobvr, n, A, n, B, n, wr, wi, b, Z2, n, Z1, n, &query, &lwork, &info);
		//compute
		lwork = int32_t(query);
		double* work = new double[lwork];
		dggev_(&jobvl, &jobvr, n, A, n, B, n, wr, wi, b, Z2, n, Z1, n, work, &lwork, &info);
		//division
		for(uint32_t i = 0; i < m_order; i++)
		{
			m_eigenvalues[0][i] /= b[i];
			m_eigenvalues[1][i] /= b[i];
		}
		//delete
		delete[] b;
		delete[] A;
		delete[] B;
		delete[] work;
		//return
		return info == 0;
	}
	bool Eigen::compute_symmetric_std_full(void)
	{
		//data
		double query;
		int32_t info;
		int32_t lwork = -1;
		const char uplo = 'U';
		const char jobz = !m_eigenvectors_computation[0] ? 'N' : 'V';
		memcpy(m_eigenvectors[0], m_data[0], m_order * m_order * sizeof(double));
		//query
		dsyev_(&jobz, &uplo, &m_order, m_eigenvectors[0], &m_order, m_eigenvalues[0], &query, &lwork, &info);
		//compute
		lwork = int32_t(query);
		double* work = new double[lwork];
		dsyev_(&jobz, &uplo, &m_order, m_eigenvectors[0], &m_order, m_eigenvalues[0], work, &lwork, &info);
		//delete
		delete[] work;
		//return
		return info == 0;
	}
	bool Eigen::compute_symmetric_gen_full(void)
	{
		//data
		double query;
		int32_t info;
		int32_t lwork = -1;
		uint32_t itype = 1;
		const char uplo = 'U';
		const char jobz = !m_eigenvectors_computation[0] ? 'N' : 'V';
		memcpy(m_eigenvectors[0], m_data[0], m_order * m_order * sizeof(double));
		//query
		double* B = new double[m_order * m_order];
		memcpy(B, m_data[1], m_order * m_order * sizeof(double));
		dsygv_(&itype, &jobz, &uplo, &m_order, m_eigenvectors[0], &m_order, B, &m_order, m_eigenvalues[0], &query, &lwork, &info);
		//compute
		lwork = int32_t(query);
		double* work = new double[lwork];
		dsygv_(&itype, &jobz, &uplo, &m_order, m_eigenvectors[0], &m_order, B, &m_order, m_eigenvalues[0], work, &lwork, &info);
		//delete
		delete[] B;
		delete[] work;
		//return
		if(info != 0) printf("info: %d\n", info);
		return info == 0;
	}
	bool Eigen::compute_symmetric_std_partial(void)
	{
		//data
		double query;
		int32_t info;
		int32_t lwork = -1;
		const char uplo = 'U';
		const char jobz = !m_eigenvectors_computation[0] ? 'N' : 'V';
		const char range = m_type == Type::Index ? 'I' : 'V';
		//setup
		int32_t* ifail = new int32_t[m_order];
		int32_t* iwork = new int32_t[5 * m_order];
		double* A = new double[m_order * m_order];
		memcpy(A, m_data[0], m_order * m_order * sizeof(double));
		//query
		uint32_t* m = &m_modes;
		const double abstol = 0;
		const uint32_t* n = &m_order;
		double* w = m_eigenvalues[0];
		double* Z = m_eigenvectors[0];
		const double* v1 = &m_value_min;
		const double* v2 = &m_value_max;
		const uint32_t* i1 = &m_index_min;
		const uint32_t* i2 = &m_index_max;
		dsyevx_(&jobz, &range, &uplo, n, A, n, v1, v2, i1, i2, &abstol, m, w, Z, n, &query, &lwork, iwork, ifail, &info);
		//compute
		lwork = int32_t(query);
		double* work = new double[lwork];
		dsyevx_(&jobz, &range, &uplo, n, A, n, v1, v2, i1, i2, &abstol, m, w, Z, n, work, &lwork, iwork, ifail, &info);
		//delete
		delete[] A;
		delete[] work;
		delete[] ifail;
		delete[] iwork;
		//return
		return info == 0;
	}
	bool Eigen::compute_symmetric_gen_partial(void)
	{
		//data
		double query;
		int32_t info;
		int32_t lwork = -1;
		const char uplo = 'U';
		const char range = m_type == Type::Index ? 'I' : 'V';
		const char jobz = !m_eigenvectors_computation[0] ? 'N' : 'V';
		//setup
		int32_t* ifail = new int32_t[m_order];
		int32_t* iwork = new int32_t[5 * m_order];
		double* A = new double[m_order * m_order];
		double* B = new double[m_order * m_order];
		memcpy(A, m_data[0], m_order * m_order * sizeof(double));
		memcpy(B, m_data[1], m_order * m_order * sizeof(double));
		//query
		uint32_t* m = &m_modes;
		const double abstol = 0;
		const uint32_t itype = 1;
		const uint32_t* n = &m_order;
		double* w = m_eigenvalues[0];
		double* Z = m_eigenvectors[0];
		const double* v1 = &m_value_min;
		const double* v2 = &m_value_max;
		const uint32_t* i1 = &m_index_min;
		const uint32_t* i2 = &m_index_max;
		dsygvx_(&itype, &jobz, &range, &uplo, n, A, n, B, n, v1, v2, i1, i2, &abstol, m, w, Z, n, &query, &lwork, iwork, ifail, &info);
		//compute
		lwork = int32_t(query);
		double* work = new double[lwork];
		dsygvx_(&itype, &jobz, &range, &uplo, n, A, n, B, n, v1, v2, i1, i2, &abstol, m, w, Z, n, work, &lwork, iwork, ifail, &info);
		//delete
		delete[] A;
		delete[] B;
		delete[] work;
		delete[] ifail;
		delete[] iwork;
		//return
		return info == 0;
	}
}