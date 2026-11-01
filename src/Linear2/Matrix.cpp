//std
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>

//Math
#include "Math/inc/Linear2/Vector.hpp"
#include "Math/inc/Linear2/Matrix.hpp"

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
	namespace linear
	{
		//constructors
		Matrix::Matrix(const Matrix& M) : Matrix{M.m_rows, M.m_cols}
		{
			memcpy(m_data_ptr, M.m_data_ref, m_rows * m_cols * sizeof(double));
		}
		Matrix::Matrix(uint32_t rows, uint32_t cols) : m_own{true}, m_rows{rows}, m_cols{cols}
		{
			m_data_ref = m_data_ptr = new double[rows * cols];
		}
		Matrix::Matrix(double* ptr, uint32_t rows, uint32_t cols) :
			m_own{false}, m_rows{rows}, m_cols{cols}, m_data_ptr{ptr}, m_data_ref{ptr}
		{
			return;
		}
		Matrix::Matrix(const double* ref, uint32_t rows, uint32_t cols) :
			m_own{false}, m_rows{rows}, m_cols{cols}, m_data_ptr{nullptr}, m_data_ref{ref}
		{
			return;
		}
		Matrix::Matrix(std::initializer_list<std::initializer_list<double>> list) : m_own{true}
		{
			//data
			const std::initializer_list<double>* data = list.begin();
			//check
			if(list.size() == 0)
			{
				throw std::runtime_error("Matrix constructor has incompatible dimensions!");
			}
			for(uint32_t i = 1; i < list.size(); i++)
			{
				if(data[i].size() != data[0].size())
				{
					throw std::runtime_error("Matrix constructor has incompatible dimensions!");
				}
			}
			//allocate
			m_rows = list.size();
			m_cols = data[0].size();
			m_data_ref = m_data_ptr = new double[m_rows * m_cols];
			//assign
			for(uint32_t i = 0; i < m_rows; i++)
			{
				const double* row_data = data[i].begin();
				for(uint32_t j = 0; j < m_cols; j++)
				{
					m_data_ptr[i + m_rows * j] = row_data[j];
				}
			}
		}
	
		//destructor
		Matrix::~Matrix(void)
		{
			if(m_own) delete[] m_data_ptr;
		}
	
		//operators
		Matrix Matrix::operator+(void) const
		{
			return *this;
		}
		Matrix Matrix::operator-(void) const
		{
			return -1 * *this;
		}
		Matrix Matrix::operator/(double s) const
		{
			return Matrix(*this) /= s;
		}
	
		Matrix Matrix::operator+(const Matrix& m) const
		{
			return Matrix(*this) += m;
		}
		Matrix Matrix::operator-(const Matrix& m) const
		{
			return Matrix(*this) -= m;
		}
		Vector Matrix::operator*(const Vector& v) const
		{
			const Matrix& m = v;
			return Vector(operator*(m));
		}
		Matrix Matrix::operator*(const Matrix& m) const
		{
			//check
			if(m_cols != m.m_rows)
			{
				throw std::runtime_error("Matrix product has incompatible dimensions!");
			}
			//compute
			Matrix r(m_rows, m.m_cols);
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
	
		Matrix& Matrix::operator=(double s)
		{
			for(uint32_t i = 0; i < m_rows * m_cols; i++)
			{
				m_data_ptr[i] = s;
			}
			return *this;
		}
		Matrix& Matrix::operator=(const Matrix& m)
		{
			//check
			if(m_rows != m.m_rows || m_cols != m.m_cols)
			{
				throw std::runtime_error("Matrix assign has incompatible dimensions!");
			}
			//assign
			memcpy(m_data_ptr, m.m_data_ref, m_rows * m_cols * sizeof(double));
			//return
			return *this;
		}
		Matrix& Matrix::operator=(const double* ref)
		{
			memcpy(m_data_ptr, ref, m_rows * m_cols * sizeof(double));
			return *this;
		}
		Matrix& Matrix::operator=(std::initializer_list<double> list)
		{
			//check
			if(list.size() != m_rows * m_cols)
			{
				throw std::runtime_error("Matrix assign has incompatible dimensions!");
			}
			//assign
			memcpy(m_data_ptr, std::data(list), m_rows * m_cols * sizeof(double));
			//return
			return *this;
		}
		Matrix& Matrix::operator=(std::initializer_list<std::initializer_list<double>> list)
		{
			//check
			if(list.size() != m_cols)
			{
				throw std::runtime_error("Matrix assign has incompatible dimensions!");
			}
			const std::initializer_list<double> *data = std::data(list);
			for(uint32_t i = 0; i < list.size(); i++)
			{
				if(data[i].size() != m_rows)
				{
					throw std::runtime_error("Matrix assign has incompatible dimensions!");
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
	
		Matrix& Matrix::operator+=(double s)
		{
			for(uint32_t i = 0; i < m_rows * m_cols; i++)
			{
				m_data_ptr[i] += s;
			}
			return *this;
		}
		Matrix& Matrix::operator-=(double s)
		{
			for(uint32_t i = 0; i < m_rows * m_cols; i++)
			{
				m_data_ptr[i] -= s;
			}
			return *this;
		}
		Matrix& Matrix::operator*=(double s)
		{
			for(uint32_t i = 0; i < m_rows * m_cols; i++)
			{
				m_data_ptr[i] *= s;
			}
			return *this;
		}
		Matrix& Matrix::operator/=(double s)
		{
			for(uint32_t i = 0; i < m_rows * m_cols; i++)
			{
				m_data_ptr[i] /= s;
			}
			return *this;
		}
	
		Matrix& Matrix::operator+=(const Matrix& m)
		{
			//check
			if(m_rows != m.m_rows || m_cols != m.m_cols)
			{
				throw std::runtime_error("Matrix increment has incompatible dimensions!");
			}
			//compute
			for(uint32_t i = 0; i < m_rows * m_cols; i++)
			{
				m_data_ptr[i] += m.m_data_ref[i];
			}
			//return
			return *this;
		}
		Matrix& Matrix::operator-=(const Matrix& m)
		{
			//check
			if(m_rows != m.m_rows || m_cols != m.m_cols)
			{
				throw std::runtime_error("Matrix increment has incompatible dimensions!");
			}
			//compute
			for(uint32_t i = 0; i < m_rows * m_cols; i++)
			{
				m_data_ptr[i] -= m.m_data_ref[i];
			}
			//return
			return *this;
		}
	
		Matrix& Matrix::operator+=(const double* ref)
		{
			for(uint32_t i = 0; i < m_rows * m_cols; i++)
			{
				m_data_ptr[i] += ref[i];
			}
			return *this;
		}
		Matrix& Matrix::operator-=(const double* ref)
		{
			for(uint32_t i = 0; i < m_rows * m_cols; i++)
			{
				m_data_ptr[i] -= ref[i];
			}
			return *this;
		}
	
		double& Matrix::operator[](uint32_t i)
		{
			return m_data_ptr[i];
		}
		double& Matrix::operator()(uint32_t i)
		{
			return m_data_ptr[i];
		}
		double& Matrix::operator()(uint32_t i, uint32_t j)
		{
			return m_data_ptr[i + m_rows * j];
		}
	
		const double& Matrix::operator[](uint32_t i) const
		{
			return m_data_ref[i];
		}
		const double& Matrix::operator()(uint32_t i) const
		{
			return m_data_ref[i];
		}
		const double& Matrix::operator()(uint32_t i, uint32_t j) const
		{
			return m_data_ref[i + m_rows * j];
		}
	
		//data
		uint32_t Matrix::rows(void) const
		{
			return m_rows;
		}
		uint32_t Matrix::cols(void) const
		{
			return m_cols;
		}
	
		double* Matrix::data(void)
		{
			return m_data_ptr;
		}
		const double* Matrix::data(void) const
		{
			return m_data_ref;
		}
	
		//bounds
		double Matrix::min(bool q, uint32_t* p) const
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
		double Matrix::max(bool q, uint32_t* p) const
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
	
		//util
		Matrix& Matrix::eye(void)
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
		Matrix& Matrix::zeros(void)
		{
			memset(m_data_ptr, 0, m_rows * m_cols * sizeof(double));
			return *this;
		}
		Matrix& Matrix::randu(double a, double b)
		{
			for(uint32_t i = 0; i < m_rows * m_cols; i++)
			{
				m_data_ptr[i] = (b - a) * rand() / RAND_MAX + a;
			}
			return *this;
		}
	
		//info
		void Matrix::print(const char* s, double v) const
		{
			if(strlen(s) != 0)
			{
				printf("%s\n", s);
			}
			for(uint32_t i = 0; i < m_rows; i++)
			{
				for(uint32_t j = 0; j < m_cols; j++)
				{
					if(v != 0 && fabs(m_data_ref[i + m_rows * j]) < v)
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
		double Matrix::norm(void) const
		{
			double s = 0;
			for(uint32_t i = 0; i < m_rows * m_cols; i++)
			{
				s += m_data_ref[i] * m_data_ref[i];
			}
			return sqrt(s);
		}
		double Matrix::trace(void) const
		{
			if(m_rows == m_cols)
			{
				double s = 0;
				for(uint32_t i = 0; i < m_rows; i++)
				{
					s += m_data_ref[i + m_rows * i];
				}
				return s;
			}
			return 0;
		}
		double Matrix::determinant(void) const
		{
			//check
			if(m_rows != m_cols)
			{
				throw std::runtime_error("Matrix determinant called on non-square Matrix!");
			}
			//data
			int32_t info;
			uint32_t* ipiv = new uint32_t[m_rows];
			double* A = new double[m_rows * m_cols];
			memcpy(A, m_data_ref, m_rows * m_cols * sizeof(double));
			//decompose
			dgetrf_(&m_rows, &m_cols, A, &m_rows, ipiv, &info);
			if(info < 0) throw std::runtime_error("Matrix determinant computation failed!");
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
		double Matrix::bilinear(const double* v) const
		{
			return bilinear(v, v);
		}
		double Matrix::bilinear(const Vector& v) const
		{
			return bilinear(v, v);
		}
		double Matrix::bilinear(const double* v1, const double* v2) const
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
		double Matrix::bilinear(const Vector& v1, const Vector& v2) const
		{
			//check
			if(m_rows != v1.m_rows || m_cols != v2.m_rows)
			{
				throw std::runtime_error("Matrix bilinear called with inconsistent dimensions!");
			}
			//return
			return bilinear(v1.m_data_ref, v2.m_data_ref);
		}
	
		Matrix Matrix::transpose(void) const
		{
			Matrix M(m_cols, m_rows);
			for(uint32_t i = 0; i < m_rows; i++)
			{
				for(uint32_t j = 0; j < m_cols; j++)
				{
					M.m_data_ptr[j + m_cols * i] = m_data_ref[i + m_rows * j];
				}
			}
			return M;
		}
		Matrix Matrix::inverse(bool* test) const
		{
			//data
			double query;
			Matrix A(*this);
			int32_t info, lwork = -1;
			//check
			if(m_rows != m_cols)
			{
				throw std::runtime_error("Matrix inverse called on a non-square Matrix!");
			}
			//query
			uint32_t* ipiv = new uint32_t[m_rows];
			dgetrf_(&m_rows, &m_cols, A.m_data_ptr, &m_rows, ipiv, &info);
			dgetri_(&m_rows, A.m_data_ptr, &m_rows, ipiv, &query, &lwork, &info);
			//inverse
			lwork = int32_t(query);
			double* work = new double[lwork];
			dgetri_(&m_rows, A.m_data_ptr, &m_rows, ipiv, work, &lwork, &info);
			//delete
			delete[] work;
			delete[] ipiv;
			if(test) *test = info == 0;
			//return
			return A;
		}
		bool Matrix::solve(Matrix& x, const Matrix& f) const
		{
			//check
			if(m_rows != m_cols)
			{
				throw std::runtime_error("Matrix solve called on a non-square Matrix!");
			}
			if(m_cols != x.m_rows || m_rows != f.m_rows || x.m_cols != f.m_cols)
			{
				throw std::runtime_error("Matrix solve called with incompatible matrices!");
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
		bool Matrix::solve(double* x, const double* f, uint32_t cols) const
		{
			//check
			if(m_rows != m_cols)
			{
				throw std::runtime_error("Matrix solve called on a non-square Matrix!");
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
	
		bool Matrix::solve_decompose(uint32_t* pivot)
		{
			//check
			if(m_rows != m_cols)
			{
				throw std::runtime_error("Matrix decompose called on a non-square Matrix!");
			}
			//decompose
			int32_t info;
			dgetrf_(&m_rows, &m_cols, m_data_ptr, &m_cols, pivot, &info);
			//return
			return info == 0;
		}
		bool Matrix::solve_substitute(const uint32_t* pivot, Matrix& x)
		{
			//check
			if(x.m_rows != m_cols)
			{
				throw std::runtime_error("Matrix substitute called with incompatible matrices!");
			}
			//return
			return solve_substitute(pivot, x.m_data_ptr, x.m_cols);
		}
		bool Matrix::solve_substitute(const uint32_t* pivot, double* x, uint32_t nhrs)
		{
			//check
			if(m_rows != m_cols)
			{
				throw std::runtime_error("Matrix substitute called on a non-square Matrix!");
			}
			//substitute
			uint32_t info;
			dgetrs_("N", &m_rows, &nhrs, m_data_ptr, &m_cols, pivot, x, &m_cols, &info);
			//return
			return info == 0;
		}
		bool Matrix::solve_substitute(const uint32_t* pivot, const Matrix& f, Matrix& x)
		{
			//check
			if(f.m_rows != m_rows || x.m_rows != m_cols || f.m_cols != x.m_cols)
			{
				throw std::runtime_error("Matrix substitute called with incompatible matrices!");
			}
			//return
			x = f.m_data_ref;
			return solve_substitute(pivot, x.m_data_ptr, x.m_cols);
		}
		bool Matrix::solve_substitute(const uint32_t* pivot, const double* f, double* x, uint32_t nrhs)
		{
			Matrix(x, m_cols, nrhs) = f;
			return solve_substitute(pivot, x, nrhs);
		}
	
		// math::Span Matrix::span(uint32_t r, uint32_t c, uint32_t n, uint32_t m)
		// {
		// 	return math::Span(*this, r, c, n, m);
		// }
	
		//stats
		double Matrix::mean(void) const
		{
			double m = 0;
			for(uint32_t i = 0; i < m_rows * m_cols; i++)
			{
				m += m_data_ref[i];
			}
			return m / m_rows / m_cols;
		}
		double Matrix::variance(void) const
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
		Matrix Matrix::eye(uint32_t rows, uint32_t cols)
		{
			return Matrix(rows, cols).eye();
		}
		Matrix Matrix::zeros(uint32_t rows, uint32_t cols)
		{
			return Matrix(rows, cols).zeros();
		}
	
		//friends
		Matrix operator*(double s, const Matrix& m)
		{
			return Matrix(m) *= s;
		}
	}
}