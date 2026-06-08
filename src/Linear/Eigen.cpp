//std
#include <cfloat>
#include <cstring>

//Math
#include "Math/inc/Linear/Eigen.hpp"

extern "C"
{
	void dsyev_(const char*, const char*, const uint32_t*, double*, const uint32_t*, double*, double*, const int32_t*, int32_t*);
	void dsyevx_(const char* jobz, const char* range, const char* uplo, const uint32_t* n, double* A, const uint32_t* lda, const double* vl, const double* vu, const uint32_t* il, const uint32_t* iu, const double* abstol, uint32_t* m, double* w, double* Z, const uint32_t* ldz, double* work, const int32_t* lwork, int32_t* iwork, int32_t* ifail, int32_t* info);
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
			}
		}
		return true;
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
		int32_t status;
		int32_t lwork = -1;
		const char uplo = 'U';
		const char jobz = !eigenvectors ? 'N' : 'V';
		memcpy(m_eigenvectors[0], m_data[0], m_order * m_order * sizeof(double));
		//query
		dsyev_(&jobz, &uplo, &m_order, m_eigenvectors[0], &m_order, m_eigenvalues[0], &query, &lwork, &status);
		//compute
		lwork = int32_t(query);
		double* work = new double[lwork];
		dsyev_(&jobz, &uplo, &m_order, m_eigenvectors[0], &m_order, m_eigenvalues[0], work, &lwork, &status);
		//delete
		delete[] work;
		//return
		return status == 0;
	}
	bool Eigen::compute_symmetric_std_partial(bool eigenvectors)
	{
		//data
		double query;
		int32_t status;
		int32_t lwork = -1;
		const double tol = 0;
		const char uplo = 'U';
		const char jobz = !eigenvectors ? 'N' : 'V';
		const char range = m_type == Type::Index ? 'I' : 'V';
		//setup
		int32_t* ifail = new int32_t[m_order];
		int32_t* iwork = new int32_t[5 * m_order];
		double* A = new double[m_order * m_order];
		memcpy(A, m_data[0], m_order * m_order * sizeof(double));
		//query
		dsyevx_(&jobz, &range, &uplo, &m_order, A, &m_order, &m_value_min, &m_value_max, &m_index_min, &m_index_max, &tol, &m_order, m_eigenvectors[0], m_eigenvectors[0], &m_order, &query, &lwork, iwork, ifail, &status);
		//delete
		delete[] A;
		delete[] ifail;
		delete[] iwork;
		//return
		return status == 0;
	}
}