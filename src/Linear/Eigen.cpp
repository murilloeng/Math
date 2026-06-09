//std
#include <cfloat>
#include <cstdio>
#include <cstring>

//Math
#include "Math/inc/Linear/Eigen.hpp"

extern "C"
{
	void dsyev_(const char*, const char*, const uint32_t*, double*, const uint32_t*, double*, double*, const int32_t*, int32_t*);
	void dsygv_(const uint32_t*, const char*, const char*, const uint32_t*, double*, const uint32_t*, double*, const uint32_t*, double*, double*, const int32_t*, int32_t*);
	void dsyevx_(const char*, const char*, const char*, const uint32_t*, double*, const uint32_t*, const double*, const double*, const uint32_t*, const uint32_t*, const double*, uint32_t*, double*, double*, const uint32_t*, double*, const int32_t*, int32_t*, int32_t*, int32_t*);
}

namespace math
{
	//constructor
	Eigen::Eigen(void) : 
		m_type{Type::Full}, m_symmetry{true}, m_order{0},
		m_value_min{-DBL_MAX}, m_value_max{+DBL_MAX}, m_index_min{0}, m_index_max{UINT32_MAX},
		m_data{nullptr, nullptr}, m_eigenvalues{nullptr, nullptr}, m_eigenvectors{nullptr, nullptr}
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

	bool Eigen::symmetry(void) const
	{
		return m_symmetry;
	}
	bool Eigen::symmetry(bool symmetry)
	{
		return m_symmetry = symmetry;
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

	//compute
	bool Eigen::compute(bool eigenvectors)
	{
		cleanup();
		allocate();
		if(m_symmetry)
		{
			if(m_data[1] == nullptr)
			{
				if(m_type == Type::Full)
				{
					return compute_symmetric_std_full(eigenvectors);
				}
				else
				{
					return compute_symmetric_std_partial(eigenvectors);
				}
			}
			else
			{
				if(m_type == Type::Full)
				{
					return compute_symmetric_gen_full(eigenvectors);
				}
			}
		}
		return true;
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
		if(!m_symmetry) m_eigenvalues[1] = new double[m_order];
		if(!m_symmetry) m_eigenvectors[1] = new double[m_order * m_order];
	}

	//compute
	bool Eigen::compute_symmetric_std_full(bool eigenvectors)
	{
		//data
		double query;
		int32_t info;
		int32_t lwork = -1;
		const char uplo = 'U';
		const char jobz = !eigenvectors ? 'N' : 'V';
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
	bool Eigen::compute_symmetric_gen_full(bool eigenvectors)
	{
		//data
		double query;
		int32_t info;
		int32_t lwork = -1;
		uint32_t itype = 1;
		const char uplo = 'U';
		const char jobz = !eigenvectors ? 'N' : 'V';
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
	bool Eigen::compute_symmetric_std_partial(bool eigenvectors)
	{
		//data
		double query;
		int32_t info;
		int32_t lwork = -1;
		const char uplo = 'U';
		const char jobz = !eigenvectors ? 'N' : 'V';
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
	bool Eigen::compute_symmetric_gen_partial(bool eigenvectors)
	{
		return true;
	}
	bool Eigen::compute_non_symmetric_std_full(bool eigenvectors)
	{
		return true;
	}
	bool Eigen::compute_non_symmetric_gen_full(bool eigenvectors)
	{
		return true;
	}
	bool Eigen::compute_non_symmetric_std_partial(bool eigenvectors)
	{
		return true;
	}
	bool Eigen::compute_non_symmetric_gen_partial(bool eigenvectors)
	{
		return true;
	}
}