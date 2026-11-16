//std
#include <cmath>
#include <cstring>
#include <algorithm>
#include <stdexcept>

//Math
#include "Math/inc/Linear/SpanX.hpp"
#include "Math/inc/Linear/MatrixX.hpp"
#include "Math/inc/Linear/VectorX.hpp"

//extern
extern "C"
{
	//inverse
	void dgetri_(const uint32_t*, double*, const uint32_t*, const uint32_t*, double*, const int32_t*, int32_t*);
	//solve
	void dgetrf_(const uint32_t*, const uint32_t*, double*, const uint32_t*, uint32_t*, int32_t*);
	void dgesv_(const uint32_t*, const uint32_t*, double*, const uint32_t*, uint32_t*, double*, const uint32_t*, int32_t*);
	void dgetrs_(const char*, const uint32_t*, const uint32_t*, double*, const uint32_t*, const uint32_t*, double*, const uint32_t*, uint32_t*);
}

namespace math
{
	//constructors
	MatrixX::MatrixX(const MatrixX& M) : MatrixX{M.m_rows, M.m_cols}
	{
		memcpy(m_data_ptr, M.m_data_ref, m_rows * m_cols * sizeof(double));
	}
	MatrixX::MatrixX(uint32_t rows, uint32_t cols) : m_own{true}, m_rows{rows}, m_cols{cols}
	{
		if(rows * cols <= MATRIX_STATIC_SIZE)
		{
			m_data_ref = m_data_ptr = m_data_mem;
		}
		else
		{
			m_data_ref = m_data_ptr = new double[rows * cols];
		}
	}
	MatrixX::MatrixX(uint32_t rows, uint32_t cols, Mode mode) : MatrixX(rows, cols)
	{
		memset(m_data_ptr, 0, m_rows * m_cols * sizeof(double));
		if(mode == Mode::eye) for(uint32_t i = 0; i < std::min(m_rows, m_cols); i++) m_data_ptr[i + m_rows * i] = 1;
	}
	MatrixX::MatrixX(double* ptr, uint32_t rows, uint32_t cols) :
		m_own{false}, m_rows{rows}, m_cols{cols}, m_data_ptr{ptr}, m_data_ref{ptr}
	{
		return;
	}
	MatrixX::MatrixX(const double* ref, uint32_t rows, uint32_t cols) :
		m_own{false}, m_rows{rows}, m_cols{cols}, m_data_ptr{nullptr}, m_data_ref{ref}
	{
		return;
	}
	MatrixX::MatrixX(std::initializer_list<std::initializer_list<double>> list) : m_own{true}
	{
		//data
		const std::initializer_list<double>* data = std::data(list);
		//check
		if(list.size() == 0)
		{
			throw std::runtime_error("MatrixX constructor has incompatible dimensions!");
		}
		for(uint32_t i = 1; i < list.size(); i++)
		{
			if(data[i].size() != data[0].size())
			{
				throw std::runtime_error("MatrixX constructor has incompatible dimensions!");
			}
		}
		//sizes
		m_cols = data[0].size();
		m_rows = list.size();
		//allocate
		if(m_rows * m_cols <= MATRIX_STATIC_SIZE)
		{
			m_data_ref = m_data_ptr = m_data_mem;
		}
		else
		{
			m_data_ref = m_data_ptr = new double[m_rows * m_cols];
		}
		//assign
		for(uint32_t i = 0; i < m_rows; i++)
		{
			const double* row_data = std::data(data[i]);
			for(uint32_t j = 0; j < m_cols; j++)
			{
				m_data_ptr[i + m_rows * j] = row_data[j];
			}
		}
	}

	//destructor
	MatrixX::~MatrixX(void)
	{
		if(m_own && m_data_ptr != m_data_mem)
		{
			delete[] m_data_ptr;
		}
	}

	//operators
	MatrixX MatrixX::operator+(void) const
	{
		return *this;
	}
	MatrixX MatrixX::operator-(void) const
	{
		return -1 * *this;
	}
	MatrixX MatrixX::operator/(double s) const
	{
		return MatrixX(*this) /= s;
	}

	MatrixX MatrixX::operator+(const MatrixX& m) const
	{
		return MatrixX(*this) += m;
	}
	MatrixX MatrixX::operator-(const MatrixX& m) const
	{
		return MatrixX(*this) -= m;
	}
	VectorX MatrixX::operator*(const VectorX& v) const
	{
		//check
		if(m_cols != v.m_rows)
		{
			throw std::runtime_error("MatrixX-VectorX product has incompatible dimensions!");
		}
		//compute
		VectorX r(m_rows);
		for(uint32_t i = 0; i < m_rows; i++)
		{
			r.m_data_ptr[i] = 0;
			for(uint32_t j = 0; j < m_cols; j++)
			{
				r.m_data_ptr[i] += m_data_ref[i + m_rows * j] * v.m_data_ref[j];
			}
		}
		return r;
	}
	MatrixX MatrixX::operator*(const MatrixX& m) const
	{
		//check
		if(m_cols != m.m_rows)
		{
			throw std::runtime_error("MatrixX product has incompatible dimensions!");
		}
		//compute
		MatrixX r(m_rows, m.m_cols);
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m.m_cols; j++)
			{
				r(i, j) = 0;
				for(uint32_t k = 0; k < m_cols; k++)
				{
					r.m_data_ptr[i + r.m_rows * j] += m_data_ref[i + m_rows * k] * m.m_data_ref[k + m.m_rows * j];
				}
			}
		}
		return r;
	}

	MatrixX& MatrixX::operator=(double s)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] = s;
		}
		return *this;
	}
	MatrixX& MatrixX::operator=(const MatrixX& m)
	{
		//check
		if(m_rows != m.m_rows || m_cols != m.m_cols)
		{
			throw std::runtime_error("MatrixX assign has incompatible dimensions!");
		}
		//assign
		memcpy(m_data_ptr, m.m_data_ref, m_rows * m_cols * sizeof(double));
		//return
		return *this;
	}
	MatrixX& MatrixX::operator=(const double* ref)
	{
		memcpy(m_data_ptr, ref, m_rows * m_cols * sizeof(double));
		return *this;
	}
	MatrixX& MatrixX::operator=(std::initializer_list<double> list)
	{
		//check
		if(list.size() != m_rows * m_cols)
		{
			throw std::runtime_error("MatrixX assign has incompatible dimensions!");
		}
		//assign
		memcpy(m_data_ptr, std::data(list), m_rows * m_cols * sizeof(double));
		//return
		return *this;
	}
	MatrixX& MatrixX::operator=(std::initializer_list<std::initializer_list<double>> list)
	{
		//check
		if(list.size() != m_cols)
		{
			throw std::runtime_error("MatrixX assign has incompatible dimensions!");
		}
		const std::initializer_list<double> *data = std::data(list);
		for(uint32_t i = 0; i < list.size(); i++)
		{
			if(data[i].size() != m_rows)
			{
				throw std::runtime_error("MatrixX assign has incompatible dimensions!");
			}
		}
		//assign
		for(uint32_t i = 0; i < m_cols; i++)
		{
			memcpy(m_data_ptr + i * m_rows, std::data(data[i]), m_rows * sizeof(double));
		}
		//return
		return *this;
	}

	MatrixX& MatrixX::operator+=(double s)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] += s;
		}
		return *this;
	}
	MatrixX& MatrixX::operator-=(double s)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] -= s;
		}
		return *this;
	}
	MatrixX& MatrixX::operator*=(double s)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] *= s;
		}
		return *this;
	}
	MatrixX& MatrixX::operator/=(double s)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] /= s;
		}
		return *this;
	}

	MatrixX& MatrixX::operator+=(const MatrixX& m)
	{
		//check
		if(m_rows != m.m_rows || m_cols != m.m_cols)
		{
			throw std::runtime_error("MatrixX increment has incompatible dimensions!");
		}
		//compute
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] += m.m_data_ref[i];
		}
		//return
		return *this;
	}
	MatrixX& MatrixX::operator-=(const MatrixX& m)
	{
		//check
		if(m_rows != m.m_rows || m_cols != m.m_cols)
		{
			throw std::runtime_error("MatrixX increment has incompatible dimensions!");
		}
		//compute
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] -= m.m_data_ref[i];
		}
		//return
		return *this;
	}

	MatrixX& MatrixX::operator+=(const double* ref)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] += ref[i];
		}
		return *this;
	}
	MatrixX& MatrixX::operator-=(const double* ref)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] -= ref[i];
		}
		return *this;
	}

	double& MatrixX::operator[](uint32_t i)
	{
		return m_data_ptr[i];
	}
	double& MatrixX::operator()(uint32_t i)
	{
		return m_data_ptr[i];
	}
	double& MatrixX::operator()(uint32_t i, uint32_t j)
	{
		return m_data_ptr[i + m_rows * j];
	}

	const double& MatrixX::operator[](uint32_t i) const
	{
		return m_data_ref[i];
	}
	const double& MatrixX::operator()(uint32_t i) const
	{
		return m_data_ref[i];
	}
	const double& MatrixX::operator()(uint32_t i, uint32_t j) const
	{
		return m_data_ref[i + m_rows * j];
	}

	//data
	uint32_t MatrixX::rows(void) const
	{
		return m_rows;
	}
	uint32_t MatrixX::cols(void) const
	{
		return m_cols;
	}

	double* MatrixX::data(void)
	{
		return m_data_ptr;
	}
	const double* MatrixX::data(void) const
	{
		return m_data_ref;
	}

	//bounds
	double MatrixX::min(bool q, uint32_t* p) const
	{
		if(p) *p = 0;
		double w, v = q ? fabs(m_data_ref[0]) : m_data_ref[0];
		for(uint32_t i = 1; i < m_rows * m_cols; i++)
		{
			w = q ? fabs(m_data_ref[i]) : m_data_ref[i];
			if(w < v)
			{
				v = w;
				if(p) *p = i;
			}
		}
		return v;
	}
	double MatrixX::max(bool q, uint32_t* p) const
	{
		if(p) *p = 0;
		double w, v = q ? fabs(m_data_ref[0]) : m_data_ref[0];
		for(uint32_t i = 1; i < m_rows * m_cols; i++)
		{
			w = q ? fabs(m_data_ref[i]) : m_data_ref[i];
			if(w > v)
			{
				v = w;
				if(p) *p = i;
			}
		}
		return v;
	}

	//utility
	MatrixX& MatrixX::eye(void)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				m_data_ptr[i + m_rows * j] = i == j;
			}
		}
		return *this;
	}
	MatrixX& MatrixX::zeros(void)
	{
		memset(m_data_ptr, 0, m_rows * m_cols * sizeof(double));
		return *this;
	}
	MatrixX& MatrixX::randu(double a, double b)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] = (b - a) * rand() / RAND_MAX + a;
		}
		return *this;
	}

	//print
	void MatrixX::print(const char* header, double cutoff) const
	{
		if(strlen(header) != 0)
		{
			printf("%s\n", header);
		}
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				if(cutoff != 0 && fabs(m_data_ref[i + m_rows * j]) < cutoff)
				{
					printf("--------- ");
				}
				else
				{
					printf("%+.2e ", m_data_ref[i + m_rows * j]);
				}
			}
			printf("\n");
		}
	}

	//linear
	double MatrixX::norm(void) const
	{
		double s = 0;
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			s += m_data_ref[i] * m_data_ref[i];
		}
		return sqrt(s);
	}
	double MatrixX::trace(void) const
	{
		double s = 0;
		for(uint32_t i = 0; i < std::min(m_rows, m_cols); i++)
		{
			s += m_data_ref[i + m_rows * i];
		}
		return s;
	}
	double MatrixX::determinant(void) const
	{
		//check
		if(m_rows != m_cols)
		{
			throw std::runtime_error("MatrixX determinant called on non-square MatrixX!");
		}
		//data
		int32_t info;
		uint32_t* ipiv = new uint32_t[m_rows];
		double* A = new double[m_rows * m_cols];
		memcpy(A, m_data_ref, m_rows * m_cols * sizeof(double));
		//decompose
		dgetrf_(&m_rows, &m_cols, A, &m_rows, ipiv, &info);
		if(info < 0) throw std::runtime_error("MatrixX determinant computation failed!");
		//determinant
		double d = 1;
		for(uint32_t i = 0; i < m_rows; i++)
		{
			d *= A[i + m_rows * i];
			if(i + 1 != ipiv[i]) d *= -1;
		}
		//delete
		delete[] A;
		delete[] ipiv;
		//return
		return d;
	}
	double MatrixX::bilinear(const double* v) const
	{
		return bilinear(v, v);
	}
	double MatrixX::bilinear(const VectorX& v) const
	{
		return bilinear(v, v);
	}
	double MatrixX::bilinear(const double* v1, const double* v2) const
	{
		//compute
		double s = 0;
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				s += v1[i] * m_data_ref[i + m_rows * j] * v2[j];
			}
		}
		//return
		return s;
	}
	double MatrixX::bilinear(const VectorX& v1, const VectorX& v2) const
	{
		//check
		if(m_rows != v1.m_rows || m_cols != v2.m_rows)
		{
			throw std::runtime_error("MatrixX bilinear called with inconsistent dimensions!");
		}
		//return
		return bilinear(v1.m_data_ref, v2.m_data_ref);
	}

	MatrixX MatrixX::transpose(void) const
	{
		MatrixX M(m_cols, m_rows);
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				M.m_data_ptr[j + m_cols * i] = m_data_ref[i + m_rows * j];
			}
		}
		return M;
	}
	MatrixX MatrixX::inverse(bool* test) const
	{
		//check
		if(m_rows != m_cols)
		{
			throw std::runtime_error("MatrixX inverse called on a non-square MatrixX!");
		}
		//query
		double query;
		MatrixX A(*this);
		int32_t info, lwork = -1;
		uint32_t* ipiv = new uint32_t[m_rows];
		dgetrf_(&m_rows, &m_cols, A.m_data_ptr, &m_rows, ipiv, &info);
		dgetri_(&m_rows, A.m_data_ptr, &m_rows, ipiv, &query, &lwork, &info);
		//inverse
		lwork = int32_t(query);
		double* work = new double[lwork];
		dgetri_(&m_rows, A.m_data_ptr, &m_rows, ipiv, work, &lwork, &info);
		//delete
		delete[] ipiv;
		delete[] work;
		if(test) *test = info == 0;
		//return
		return A;
	}
	bool MatrixX::solve(MatrixX& x, const MatrixX& f) const
	{
		//check
		if(m_rows != m_cols)
		{
			throw std::runtime_error("MatrixX solve called on a non-square MatrixX!");
		}
		if(m_cols != x.m_rows || m_rows != f.m_rows || x.m_cols != f.m_cols)
		{
			throw std::runtime_error("MatrixX solve called with incompatible matrices!");
		}
		//solve
		int32_t info;
		uint32_t* ipiv = new uint32_t[m_rows];
		double* A = new double[m_rows * m_cols];
		memcpy(A, m_data_ref, m_rows * m_cols * sizeof(double));
		memcpy(x.m_data_ptr, f.m_data_ref, x.m_rows * x.m_cols * sizeof(double));
		dgesv_(&m_rows, &x.m_cols, A, &m_rows, ipiv, x.m_data_ptr, &m_rows, &info);
		//delete
		delete[] A;
		delete[] ipiv;
		//return
		return info == 0;
	}
	bool MatrixX::solve(VectorX& x, const VectorX& f) const
	{
		//check
		if(m_rows != m_cols)
		{
			throw std::runtime_error("MatrixX solve called on a non-square MatrixX!");
		}
		if(m_cols != x.m_rows || m_rows != f.m_rows)
		{
			throw std::runtime_error("MatrixX solve called with incompatible matrices!");
		}
		//solve
		int32_t info;
		const uint32_t cols = 1;
		uint32_t* ipiv = new uint32_t[m_rows];
		double* A = new double[m_rows * m_cols];
		memcpy(A, m_data_ref, m_rows * m_cols * sizeof(double));
		memcpy(x.m_data_ptr, f.m_data_ref, x.m_rows * sizeof(double));
		dgesv_(&m_rows, &cols, A, &m_rows, ipiv, x.m_data_ptr, &m_rows, &info);
		//delete
		delete[] A;
		delete[] ipiv;
		//return
		return info == 0;
	}
	bool MatrixX::solve(double* x, const double* f, uint32_t cols) const
	{
		//check
		if(m_rows != m_cols)
		{
			throw std::runtime_error("MatrixX solve called on a non-square MatrixX!");
		}
		//solve
		int32_t info;
		uint32_t* ipiv = new uint32_t[m_rows];
		double* A = new double[m_rows * m_cols];
		memcpy(x, f, m_rows * cols * sizeof(double));
		memcpy(A, m_data_ref, m_rows * m_cols * sizeof(double));
		dgesv_(&m_rows, &cols, A, &m_rows, ipiv, x, &m_rows, &info);
		//delete
		delete[] A;
		delete[] ipiv;
		//return
		return info == 0;
	}

	bool MatrixX::solve_decompose(uint32_t* pivot)
	{
		//check
		if(m_rows != m_cols)
		{
			throw std::runtime_error("MatrixX decompose called on a non-square MatrixX!");
		}
		//decompose
		int32_t info;
		dgetrf_(&m_rows, &m_cols, m_data_ptr, &m_cols, pivot, &info);
		//return
		return info == 0;
	}
	bool MatrixX::solve_substitute(const uint32_t* pivot, MatrixX& x)
	{
		//check
		if(x.m_rows != m_cols)
		{
			throw std::runtime_error("MatrixX substitute called with incompatible matrices!");
		}
		//return
		return solve_substitute(pivot, x.m_data_ptr, x.m_cols);
	}
	bool MatrixX::solve_substitute(const uint32_t* pivot, double* x, uint32_t nhrs)
	{
		//check
		if(m_rows != m_cols)
		{
			throw std::runtime_error("MatrixX substitute called on a non-square MatrixX!");
		}
		//substitute
		uint32_t info;
		dgetrs_("N", &m_rows, &nhrs, m_data_ptr, &m_cols, pivot, x, &m_cols, &info);
		//return
		return info == 0;
	}
	bool MatrixX::solve_substitute(const uint32_t* pivot, const MatrixX& f, MatrixX& x)
	{
		//check
		if(f.m_rows != m_rows || x.m_rows != m_cols || f.m_cols != x.m_cols)
		{
			throw std::runtime_error("MatrixX substitute called with incompatible matrices!");
		}
		//return
		x = f.m_data_ref;
		return solve_substitute(pivot, x.m_data_ptr, x.m_cols);
	}
	bool MatrixX::solve_substitute(const uint32_t* pivot, const double* f, double* x, uint32_t nrhs)
	{
		MatrixX(x, m_cols, nrhs) = f;
		return solve_substitute(pivot, x, nrhs);
	}

	SpanX MatrixX::span(uint32_t row, uint32_t col, uint32_t rows, uint32_t cols)
	{
		return SpanX(*this, row, col, rows, cols);
	}

	//stats
	double MatrixX::mean(void) const
	{
		double m = 0;
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m += m_data_ref[i];
		}
		return m / m_rows / m_cols;
	}
	double MatrixX::variance(void) const
	{
		double s = 0;
		const double m = mean();
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			s += m_data_ref[i] * m_data_ref[i];
		}
		return s / m_rows / m_cols - m * m;
	}

	//static
	MatrixX MatrixX::eye(uint32_t rows, uint32_t cols)
	{
		return MatrixX(rows, cols).eye();
	}
	MatrixX MatrixX::zeros(uint32_t rows, uint32_t cols)
	{
		return MatrixX(rows, cols).zeros();
	}

	//friends
	MatrixX operator*(double s, const MatrixX& m)
	{
		return MatrixX(m) *= s;
	}
}