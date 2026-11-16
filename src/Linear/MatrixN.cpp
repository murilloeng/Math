//std
#include <cmath>
#include <cstring>
#include <algorithm>
#include <stdexcept>

//Math
#include "Math/inc/Linear/SpanN.hpp"
#include "Math/inc/Linear/MatrixN.hpp"
#include "Math/inc/Linear/VectorN.hpp"

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
	MatrixN::MatrixN(const MatrixN& M) : MatrixN{M.m_rows, M.m_cols}
	{
		memcpy(m_data_ptr, M.m_data_ref, m_rows * m_cols * sizeof(double));
	}
	MatrixN::MatrixN(uint32_t rows, uint32_t cols) : m_own{true}, m_rows{rows}, m_cols{cols}
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
	MatrixN::MatrixN(uint32_t rows, uint32_t cols, Mode mode) : MatrixN(rows, cols)
	{
		memset(m_data_ptr, 0, m_rows * m_cols * sizeof(double));
		if(mode == Mode::eye) for(uint32_t i = 0; i < std::min(m_rows, m_cols); i++) m_data_ptr[i + m_rows * i] = 1;
	}
	MatrixN::MatrixN(double* ptr, uint32_t rows, uint32_t cols) :
		m_own{false}, m_rows{rows}, m_cols{cols}, m_data_ptr{ptr}, m_data_ref{ptr}
	{
		return;
	}
	MatrixN::MatrixN(const double* ref, uint32_t rows, uint32_t cols) :
		m_own{false}, m_rows{rows}, m_cols{cols}, m_data_ptr{nullptr}, m_data_ref{ref}
	{
		return;
	}
	MatrixN::MatrixN(std::initializer_list<std::initializer_list<double>> list) : m_own{true}
	{
		//data
		const std::initializer_list<double>* data = std::data(list);
		//check
		if(list.size() == 0)
		{
			throw std::runtime_error("MatrixN constructor has incompatible dimensions!");
		}
		for(uint32_t i = 1; i < list.size(); i++)
		{
			if(data[i].size() != data[0].size())
			{
				throw std::runtime_error("MatrixN constructor has incompatible dimensions!");
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
	MatrixN::~MatrixN(void)
	{
		if(m_own && m_data_ptr != m_data_mem)
		{
			delete[] m_data_ptr;
		}
	}

	//operators
	MatrixN MatrixN::operator+(void) const
	{
		return *this;
	}
	MatrixN MatrixN::operator-(void) const
	{
		return -1 * *this;
	}
	MatrixN MatrixN::operator/(double s) const
	{
		return MatrixN(*this) /= s;
	}

	MatrixN MatrixN::operator+(const MatrixN& m) const
	{
		return MatrixN(*this) += m;
	}
	MatrixN MatrixN::operator-(const MatrixN& m) const
	{
		return MatrixN(*this) -= m;
	}
	VectorN MatrixN::operator*(const VectorN& v) const
	{
		//check
		if(m_cols != v.m_rows)
		{
			throw std::runtime_error("MatrixN-VectorN product has incompatible dimensions!");
		}
		//compute
		VectorN r(m_rows);
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
	MatrixN MatrixN::operator*(const MatrixN& m) const
	{
		//check
		if(m_cols != m.m_rows)
		{
			throw std::runtime_error("MatrixN product has incompatible dimensions!");
		}
		//compute
		MatrixN r(m_rows, m.m_cols);
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

	MatrixN& MatrixN::operator=(double s)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] = s;
		}
		return *this;
	}
	MatrixN& MatrixN::operator=(const MatrixN& m)
	{
		//check
		if(m_rows != m.m_rows || m_cols != m.m_cols)
		{
			throw std::runtime_error("MatrixN assign has incompatible dimensions!");
		}
		//assign
		memcpy(m_data_ptr, m.m_data_ref, m_rows * m_cols * sizeof(double));
		//return
		return *this;
	}
	MatrixN& MatrixN::operator=(const double* ref)
	{
		memcpy(m_data_ptr, ref, m_rows * m_cols * sizeof(double));
		return *this;
	}
	MatrixN& MatrixN::operator=(std::initializer_list<double> list)
	{
		//check
		if(list.size() != m_rows * m_cols)
		{
			throw std::runtime_error("MatrixN assign has incompatible dimensions!");
		}
		//assign
		memcpy(m_data_ptr, std::data(list), m_rows * m_cols * sizeof(double));
		//return
		return *this;
	}
	MatrixN& MatrixN::operator=(std::initializer_list<std::initializer_list<double>> list)
	{
		//check
		if(list.size() != m_cols)
		{
			throw std::runtime_error("MatrixN assign has incompatible dimensions!");
		}
		const std::initializer_list<double> *data = std::data(list);
		for(uint32_t i = 0; i < list.size(); i++)
		{
			if(data[i].size() != m_rows)
			{
				throw std::runtime_error("MatrixN assign has incompatible dimensions!");
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
	template<uint32_t r, uint32_t c> MatrixN& MatrixN::operator=(const Matrix<r, c>& M)
	{
		//check
		if(m_rows != r || m_cols != c)
		{
			throw std::runtime_error("MatrixN assign has incompatible dimensions!");
		}
		//assign
		memcpy(m_data_ptr, M.data(), r * c * sizeof(double));
		//return
		return *this;
	}

	MatrixN& MatrixN::operator+=(double s)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] += s;
		}
		return *this;
	}
	MatrixN& MatrixN::operator-=(double s)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] -= s;
		}
		return *this;
	}
	MatrixN& MatrixN::operator*=(double s)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] *= s;
		}
		return *this;
	}
	MatrixN& MatrixN::operator/=(double s)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] /= s;
		}
		return *this;
	}

	MatrixN& MatrixN::operator+=(const MatrixN& m)
	{
		//check
		if(m_rows != m.m_rows || m_cols != m.m_cols)
		{
			throw std::runtime_error("MatrixN increment has incompatible dimensions!");
		}
		//compute
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] += m.m_data_ref[i];
		}
		//return
		return *this;
	}
	MatrixN& MatrixN::operator-=(const MatrixN& m)
	{
		//check
		if(m_rows != m.m_rows || m_cols != m.m_cols)
		{
			throw std::runtime_error("MatrixN increment has incompatible dimensions!");
		}
		//compute
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] -= m.m_data_ref[i];
		}
		//return
		return *this;
	}

	MatrixN& MatrixN::operator+=(const double* ref)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] += ref[i];
		}
		return *this;
	}
	MatrixN& MatrixN::operator-=(const double* ref)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] -= ref[i];
		}
		return *this;
	}

	double& MatrixN::operator[](uint32_t i)
	{
		return m_data_ptr[i];
	}
	double& MatrixN::operator()(uint32_t i)
	{
		return m_data_ptr[i];
	}
	double& MatrixN::operator()(uint32_t i, uint32_t j)
	{
		return m_data_ptr[i + m_rows * j];
	}

	const double& MatrixN::operator[](uint32_t i) const
	{
		return m_data_ref[i];
	}
	const double& MatrixN::operator()(uint32_t i) const
	{
		return m_data_ref[i];
	}
	const double& MatrixN::operator()(uint32_t i, uint32_t j) const
	{
		return m_data_ref[i + m_rows * j];
	}

	//data
	uint32_t MatrixN::rows(void) const
	{
		return m_rows;
	}
	uint32_t MatrixN::cols(void) const
	{
		return m_cols;
	}

	double* MatrixN::data(void)
	{
		return m_data_ptr;
	}
	const double* MatrixN::data(void) const
	{
		return m_data_ref;
	}

	//bounds
	double MatrixN::min(bool q, uint32_t* p) const
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
	double MatrixN::max(bool q, uint32_t* p) const
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
	MatrixN& MatrixN::eye(void)
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
	MatrixN& MatrixN::zeros(void)
	{
		memset(m_data_ptr, 0, m_rows * m_cols * sizeof(double));
		return *this;
	}
	MatrixN& MatrixN::randu(double a, double b)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_data_ptr[i] = (b - a) * rand() / RAND_MAX + a;
		}
		return *this;
	}

	//print
	void MatrixN::print(const char* header, double cutoff) const
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
	double MatrixN::norm(void) const
	{
		double s = 0;
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			s += m_data_ref[i] * m_data_ref[i];
		}
		return sqrt(s);
	}
	double MatrixN::trace(void) const
	{
		double s = 0;
		for(uint32_t i = 0; i < std::min(m_rows, m_cols); i++)
		{
			s += m_data_ref[i + m_rows * i];
		}
		return s;
	}
	double MatrixN::determinant(void) const
	{
		//check
		if(m_rows != m_cols)
		{
			throw std::runtime_error("MatrixN determinant called on non-square MatrixN!");
		}
		//data
		int32_t info;
		uint32_t* ipiv = new uint32_t[m_rows];
		double* A = new double[m_rows * m_cols];
		memcpy(A, m_data_ref, m_rows * m_cols * sizeof(double));
		//decompose
		dgetrf_(&m_rows, &m_cols, A, &m_rows, ipiv, &info);
		if(info < 0) throw std::runtime_error("MatrixN determinant computation failed!");
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
	double MatrixN::bilinear(const double* v) const
	{
		return bilinear(v, v);
	}
	double MatrixN::bilinear(const VectorN& v) const
	{
		return bilinear(v, v);
	}
	double MatrixN::bilinear(const double* v1, const double* v2) const
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
	double MatrixN::bilinear(const VectorN& v1, const VectorN& v2) const
	{
		//check
		if(m_rows != v1.m_rows || m_cols != v2.m_rows)
		{
			throw std::runtime_error("MatrixN bilinear called with inconsistent dimensions!");
		}
		//return
		return bilinear(v1.m_data_ref, v2.m_data_ref);
	}

	MatrixN MatrixN::transpose(void) const
	{
		MatrixN M(m_cols, m_rows);
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				M.m_data_ptr[j + m_cols * i] = m_data_ref[i + m_rows * j];
			}
		}
		return M;
	}
	MatrixN MatrixN::inverse(bool* test) const
	{
		//check
		if(m_rows != m_cols)
		{
			throw std::runtime_error("MatrixN inverse called on a non-square MatrixN!");
		}
		//query
		double query;
		MatrixN A(*this);
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
	bool MatrixN::solve(MatrixN& x, const MatrixN& f) const
	{
		//check
		if(m_rows != m_cols)
		{
			throw std::runtime_error("MatrixN solve called on a non-square MatrixN!");
		}
		if(m_cols != x.m_rows || m_rows != f.m_rows || x.m_cols != f.m_cols)
		{
			throw std::runtime_error("MatrixN solve called with incompatible matrices!");
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
	bool MatrixN::solve(VectorN& x, const VectorN& f) const
	{
		//check
		if(m_rows != m_cols)
		{
			throw std::runtime_error("MatrixN solve called on a non-square MatrixN!");
		}
		if(m_cols != x.m_rows || m_rows != f.m_rows)
		{
			throw std::runtime_error("MatrixN solve called with incompatible matrices!");
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
	bool MatrixN::solve(double* x, const double* f, uint32_t cols) const
	{
		//check
		if(m_rows != m_cols)
		{
			throw std::runtime_error("MatrixN solve called on a non-square MatrixN!");
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

	bool MatrixN::solve_decompose(uint32_t* pivot)
	{
		//check
		if(m_rows != m_cols)
		{
			throw std::runtime_error("MatrixN decompose called on a non-square MatrixN!");
		}
		//decompose
		int32_t info;
		dgetrf_(&m_rows, &m_cols, m_data_ptr, &m_cols, pivot, &info);
		//return
		return info == 0;
	}
	bool MatrixN::solve_substitute(const uint32_t* pivot, MatrixN& x)
	{
		//check
		if(x.m_rows != m_cols)
		{
			throw std::runtime_error("MatrixN substitute called with incompatible matrices!");
		}
		//return
		return solve_substitute(pivot, x.m_data_ptr, x.m_cols);
	}
	bool MatrixN::solve_substitute(const uint32_t* pivot, double* x, uint32_t nhrs)
	{
		//check
		if(m_rows != m_cols)
		{
			throw std::runtime_error("MatrixN substitute called on a non-square MatrixN!");
		}
		//substitute
		uint32_t info;
		dgetrs_("N", &m_rows, &nhrs, m_data_ptr, &m_cols, pivot, x, &m_cols, &info);
		//return
		return info == 0;
	}
	bool MatrixN::solve_substitute(const uint32_t* pivot, const MatrixN& f, MatrixN& x)
	{
		//check
		if(f.m_rows != m_rows || x.m_rows != m_cols || f.m_cols != x.m_cols)
		{
			throw std::runtime_error("MatrixN substitute called with incompatible matrices!");
		}
		//return
		x = f.m_data_ref;
		return solve_substitute(pivot, x.m_data_ptr, x.m_cols);
	}
	bool MatrixN::solve_substitute(const uint32_t* pivot, const double* f, double* x, uint32_t nrhs)
	{
		MatrixN(x, m_cols, nrhs) = f;
		return solve_substitute(pivot, x, nrhs);
	}

	SpanN MatrixN::span(uint32_t row, uint32_t col, uint32_t rows, uint32_t cols)
	{
		return SpanN(*this, row, col, rows, cols);
	}

	//stats
	double MatrixN::mean(void) const
	{
		double m = 0;
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m += m_data_ref[i];
		}
		return m / m_rows / m_cols;
	}
	double MatrixN::variance(void) const
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
	MatrixN MatrixN::eye(uint32_t rows, uint32_t cols)
	{
		return MatrixN(rows, cols).eye();
	}
	MatrixN MatrixN::zeros(uint32_t rows, uint32_t cols)
	{
		return MatrixN(rows, cols).zeros();
	}

	//friends
	MatrixN operator*(double s, const MatrixN& m)
	{
		return MatrixN(m) *= s;
	}
}