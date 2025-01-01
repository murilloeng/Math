//std
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>

//math
#include "Math/Math/inc/misc/util.hpp"
#include "Math/Math/inc/linear/mat3.hpp"
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/linear/matrix.hpp"

//extern
extern "C"
{
	//inverse
	void dgetrf_(const uint32_t*, const uint32_t*, double*, const uint32_t*, uint32_t*, int32_t*);
	void dgetri_(const uint32_t*, double*, const uint32_t*, const uint32_t*, double*, const int32_t*, int32_t*);
	//solve
	void dgesv_(const uint32_t*, const uint32_t*, double*, const uint32_t*, uint32_t*, double*, const uint32_t*, int32_t*);
	//svd
	void dgesdd_(const char*, const uint32_t*, const uint32_t*, double*, const uint32_t*, double*, double*, const uint32_t*, double*, const uint32_t*, double*, const uint32_t*, int32_t*, int32_t*);
	//eigen
	void dsyev_(const char*, const char*, const uint32_t*, double*, const uint32_t*, double*, double*, int32_t*, int32_t*);
	void dsygv_(const uint32_t*, const char*, const char*, const uint32_t*, double*, const uint32_t*, double*, const uint32_t*, double*, double*, int32_t*, int32_t*);
}

namespace math
{
	//constructors
	matrix::matrix(void) : matrix(0, 0)
	{
		return;
	}
	matrix::matrix(const matrix& m) : matrix(m.m_rows, m.m_cols)
	{
		memcpy(m_ptr, m.m_ref, m_rows * m_cols * sizeof(double));
	}
	matrix::matrix(uint32_t rows, uint32_t cols, mode init) : m_own(true), m_rows(rows), m_cols(cols)
	{
		//data
		if(rows * cols <= MATRIX_STATIC_SIZE)
		{
			m_ref = m_ptr = m_mem;
		}
		else
		{
			m_ref = m_ptr = new double[rows * cols];
		}
		//setup
		if(init == mode::eye)
		{
			eye();
		}
		else if(init == mode::zeros)
		{
			zeros();
		}
	}
	matrix::matrix(double* ptr, uint32_t rows, uint32_t cols, mode init) : 
		m_own(false), m_ptr(ptr), m_rows(rows), m_cols(cols), m_ref(ptr)
	{
		if(init == mode::eye)
		{
			eye();
		}
		else if(init == mode::zeros)
		{
			zeros();
		}
	}
	matrix::matrix(const double* ref, uint32_t rows, uint32_t cols) : 
		m_own(false), m_ptr(nullptr), m_rows(rows), m_cols(cols), m_ref(ref)
	{
		return;
	}
	matrix::matrix(std::initializer_list<std::initializer_list<double>> list, bool columns) : 
		m_own(true)
	{
		//data
		const std::initializer_list<double>* data = std::data(list);
		//check
		for(uint32_t i = 0; i < list.size(); i++)
		{
			if(data[i].size() != data[0].size())
			{
				fprintf(stderr, "\tError: Matrix constructor with incompatible dimensions!\n");
				exit(EXIT_FAILURE);
			}
		}
		//sizes
		m_cols = columns ? list.size() : data[0].size();
		m_rows = columns ? data[0].size() : list.size();
		//allocate
		if(m_rows * m_cols <= MATRIX_STATIC_SIZE)
		{
			m_ref = m_ptr = m_mem;
		}
		else
		{
			m_ref = m_ptr = new double[m_rows * m_cols];
		}
		//assign
		if(columns)
		{
			for(uint32_t i = 0; i < m_cols; i++)
			{
				memcpy(m_ptr + i * m_rows, std::data(data[i]), m_rows * sizeof(double));
			}
		}
		else
		{
			for(uint32_t i = 0; i < m_rows; i++)
			{
				const double* row_data = std::data(data[i]);
				for(uint32_t j = 0; j < m_cols; j++)
				{
					m_ptr[i + m_rows * j] = row_data[j];
				}
			}
		}
	}

	//destructor
	matrix::~matrix(void)
	{
		if(m_own && m_ptr != m_mem)
		{
			delete[] m_ptr;
		}
	}

	//serialization
	void matrix::load(const char* path)
	{
		//file
		FILE* file = fopen(path, "r");
		//sizes
		if(fscanf(file, "%d %d", &m_rows, &m_cols) != 2)
		{
			printf("\tError: Unable to load matrix!\n");
			exit(EXIT_FAILURE);
		}
		//allocate
		if(m_own && m_ptr != m_mem)
		{
			delete[] m_ptr;
		}
		m_ref = m_ptr = new double[m_rows * m_cols];
		//read
		rewind(file);
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				if(fscanf(file, "%lf", &m_ptr[i + m_rows * j]) != 1)
				{
					printf("\tError: Unable to load matrix!\n");
					exit(EXIT_FAILURE);
				}
			}
		}
		//close
		fclose(file);
	}
	void matrix::save(const char* path) const
	{
		//file
		FILE* file = fopen(path, "w");
		//size
		fprintf(file, "%04d %04d\n", m_rows, m_cols);
		//save
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				fprintf(file, "%+.6e ", m_ref[i + m_rows * j]);
			}
			fprintf(file, "\n");
		}
		//close
		fclose(file);
	}

	//operators
	matrix matrix::operator+(void) const
	{
		return *this;
	}
	matrix matrix::operator-(void) const
	{
		return -1 * *this;
	}
	matrix matrix::operator/(double s) const
	{
		return matrix(*this) /= s;
	}

	matrix matrix::operator+(const matrix& m) const
	{
		return matrix(*this) += m;
	}
	matrix matrix::operator-(const matrix& m) const
	{
		return matrix(*this) -= m;
	}
	vector matrix::operator*(const vector& v) const
	{
		const matrix& m = v;
		return vector(operator*(m));
	}
	matrix matrix::operator*(const matrix& m) const
	{
		//check
		if(m_cols != m.m_rows)
		{
			fprintf(stderr, "Error: Matrix product with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		//compute
		matrix r(m_rows, m.m_cols);
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m.m_cols; j++)
			{
				r(i, j) = 0;
				for(uint32_t k = 0; k < m_cols; k++)
				{
					r.m_ptr[i + r.m_rows * j] += m_ref[i + m_rows * k] * m.m_ref[k + m.m_rows * j];
				}
			}
		}
		return r;
	}

	matrix& matrix::operator=(double s)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_ptr[i] = s;
		}
		return *this;
	}
	matrix& matrix::operator=(const matrix& m)
	{
		//check
		if(m_rows != m.m_rows || m_cols != m.m_cols)
		{
			fprintf(stderr, "Error: Matrix assign with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		//assign
		memcpy(m_ptr, m.m_ref, m_rows * m_cols * sizeof(double));
		//return
		return *this;
	}
	matrix& matrix::operator=(const double* ref)
	{
		memcpy(m_ptr, ref, m_rows * m_cols * sizeof(double));
		return *this;
	}
	matrix& matrix::operator=(std::initializer_list<double> list)
	{
		//check
		if(list.size() != m_rows * m_cols)
		{
			fprintf(stderr, "Error: Matrix assign with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		//assign
		memcpy(m_ptr, std::data(list), m_rows * m_cols * sizeof(double));
		//return
		return *this;
	}
	matrix& matrix::operator=(std::initializer_list<std::initializer_list<double>> list)
	{
		//check
		if(list.size() != m_cols)
		{
			fprintf(stderr, "Error: Matrix assign with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		const std::initializer_list<double> *data = std::data(list);
		for(uint32_t i = 0; i < list.size(); i++)
		{
			if(data[i].size() != m_rows)
			{
				fprintf(stderr, "Error: Matrix assign with incompatible dimensions!\n");
				exit(EXIT_FAILURE);
			}
		}
		//assign
		for(uint32_t i = 0; i < m_cols; i++)
		{
			memcpy(m_ptr + i * m_rows, std::data(data[i]), m_rows * sizeof(double));
		}
		//return
		return *this;
	}

	matrix& matrix::operator+=(double s)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_ptr[i] += s;
		}
		return *this;
	}
	matrix& matrix::operator-=(double s)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_ptr[i] -= s;
		}
		return *this;
	}
	matrix& matrix::operator*=(double s)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_ptr[i] *= s;
		}
		return *this;
	}
	matrix& matrix::operator/=(double s)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_ptr[i] /= s;
		}
		return *this;
	}

	matrix& matrix::operator+=(const matrix& m)
	{
		//check
		if(m_rows != m.m_rows || m_cols != m.m_cols)
		{
			fprintf(stderr, "Error: Matrix increment with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		//compute
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_ptr[i] += m.m_ref[i];
		}
		//return
		return *this;
	}
	matrix& matrix::operator-=(const matrix& m)
	{
		//check
		if(m_rows != m.m_rows || m_cols != m.m_cols)
		{
			fprintf(stderr, "Error: Matrix increment with incompatible dimensions!\n");
			exit(EXIT_FAILURE);
		}
		//compute
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_ptr[i] -= m.m_ref[i];
		}
		//return
		return *this;
	}

	matrix& matrix::operator+=(const double* ref)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_ptr[i] += ref[i];
		}
		return *this;
	}
	matrix& matrix::operator-=(const double* ref)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_ptr[i] -= ref[i];
		}
		return *this;
	}

	double& matrix::operator[](uint32_t i)
	{
		return m_ptr[i];
	}
	double& matrix::operator()(uint32_t i)
	{
		return m_ptr[i];
	}
	double& matrix::operator()(uint32_t i, uint32_t j)
	{
		return m_ptr[i + m_rows * j];
	}

	const double& matrix::operator[](uint32_t i) const
	{
		return m_ref[i];
	}
	const double& matrix::operator()(uint32_t i) const
	{
		return m_ref[i];
	}
	const double& matrix::operator()(uint32_t i, uint32_t j) const
	{
		return m_ref[i + m_rows * j];
	}

	//data
	uint32_t matrix::rows(void) const
	{
		return m_rows;
	}
	uint32_t matrix::cols(void) const
	{
		return m_cols;
	}

	double* matrix::data(void)
	{
		return m_ptr;
	}
	const double* matrix::data(void) const
	{
		return m_ref;
	}

	//size
	matrix& matrix::resize(uint32_t rows, uint32_t cols)
	{
		//check
		if(!m_own)
		{
			fprintf(stderr, "Error: Resize called on matrix that does not own memory!");
			exit(EXIT_FAILURE);
		}
		//memory
		if(rows * cols != m_rows * m_cols)
		{
			if(m_ptr != m_mem)
			{
				delete[] m_ptr;
			}
			if(rows * cols <= MATRIX_STATIC_SIZE)
			{
				m_ref = m_ptr = m_mem;
			}
			else
			{
				m_ref = m_ptr = new double[rows * cols];
			}
		}
		//setup
		m_rows = rows;
		m_cols = cols;
		//return
		return *this;
	}

	//bounds
	double matrix::min(bool q, uint32_t* p) const
	{
		if(p) *p = 0;
		double w, v = q ? fabs(m_ref[0]) : m_ref[0];
		for(uint32_t i = 1; i < m_rows * m_cols; i++)
		{
			w = q ? fabs(m_ref[i]) : m_ref[i];
			if(w < v)
			{
				v = w;
				if(p) *p = i;
			}
		}
		return v;
	}
	double matrix::max(bool q, uint32_t* p) const
	{
		if(p) *p = 0;
		double w, v = q ? fabs(m_ref[0]) : m_ref[0];
		for(uint32_t i = 1; i < m_rows * m_cols; i++)
		{
			w = q ? fabs(m_ref[i]) : m_ref[i];
			if(w > v)
			{
				v = w;
				if(p) *p = i;
			}
		}
		return v;
	}

	//util
	matrix& matrix::eye(void)
	{
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				m_ptr[i + m_rows * j] = i == j;
			}
		}
		return *this;
	}
	matrix& matrix::zeros(void)
	{
		memset(m_ptr, 0, m_rows * m_cols * sizeof(double));
		return *this;
	}
	matrix& matrix::swap_rows(uint32_t a, uint32_t b)
	{
		//check
		if(a >= m_rows || b >= m_rows)
		{
			fprintf(stderr, "Matrix: swap_rows called with index out of range!");
			exit(EXIT_FAILURE);
		}
		//swap
		if(a != b)
		{
			for(uint32_t i = 0; i < m_cols; i++)
			{
				swap(m_ptr[a + m_rows * i], m_ptr[b + m_rows * i]);
			}
		}
		return *this;
	}
	matrix& matrix::swap_cols(uint32_t a, uint32_t b)
	{
		//check
		if(a >= m_cols || b >= m_cols)
		{
			fprintf(stderr, "Matrix: swap_cols called with index out of range!");
			exit(EXIT_FAILURE);
		}
		//swap
		if(a != b)
		{
			for(uint32_t i = 0; i < m_rows; i++)
			{
				swap(m_ptr[i + m_rows * a], m_ptr[i + m_rows * b]);
			}
		}
		return *this;
	}
	matrix& matrix::randu(double a, double b)
	{
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m_ptr[i] = (b - a) * rand() / RAND_MAX + a;
		}
		return *this;
	}

	//info
	void matrix::print(const char* s, double v) const
	{
		if(strlen(s) != 0)
		{
			printf("%s\n", s);
		}
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				if(v != 0 && fabs(m_ref[i + m_rows * j]) < v)
				{
					printf("--------- ");
				}
				else
				{
					printf("%+.2e ", m_ref[i + m_rows * j]);
				}
			}
			printf("\n");
		}
	}

	//linear
	double matrix::norm(void) const
	{
		double s = 0;
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			s += m_ref[i] * m_ref[i];
		}
		return sqrt(s);
	}
	double matrix::trace(void) const
	{
		if(m_rows == m_cols)
		{
			double s = 0;
			for(uint32_t i = 0; i < m_rows; i++)
			{
				for(uint32_t j = 0; j < m_cols; j++)
				{
					s += m_ref[i + m_rows * j];
				}
			}
			return s;
		}
		return 0;
	}
	double matrix::determinant(void) const
	{
		//data
		int32_t status;
		matrix M(*this);
		uint32_t* pivot = (uint32_t*) alloca(m_rows * sizeof(uint32_t));
		//check
		if(m_rows != m_cols)
		{
			fprintf(stderr, "Error: Determinant called on non-square matrix!\n");
			exit(EXIT_FAILURE);
		}
		//decompose
		dgetrf_(&m_rows, &m_cols, M.m_ptr, &m_rows, pivot, &status);
		//check
		if(status != 0)
		{
			fprintf(stderr, "Error: Determinant failed!\n");
			exit(EXIT_FAILURE);
		}
		//determinant
		double d = 1;
		for(uint32_t i = 0; i < m_rows; i++)
		{
			d *= M.m_ptr[i + m_rows * i];
			if(i + 1 != pivot[i]) d *= -1;
		}
		//return
		return d;
	}

	matrix matrix::inverse(bool* test) const
	{
		//data
		matrix M(*this);
		double query, *work;
		int32_t status, lwork = -1;
		uint32_t* pivot = (uint32_t*) alloca(m_rows * sizeof(uint32_t));
		//check
		if(m_rows != m_cols)
		{
			if(test) *test = false;
			fprintf(stderr, "Error: Inverse called on non-square system!\n");
			exit(EXIT_FAILURE);
		}
		//query
		dgetrf_(&m_rows, &m_cols, M.m_ptr, &m_rows, pivot, &status);
		dgetri_(&m_rows, M.m_ptr, &m_rows, pivot, &query, &lwork, &status);
		//inverse
		lwork = int32_t(query);
		work = (double*) alloca(lwork * sizeof(double));
		dgetri_(&m_rows, M.m_ptr, &m_rows, pivot, work, &lwork, &status);
		//test
		if(test) *test = status == 0;
		//return
		return M;
	}
	matrix matrix::transpose(void) const
	{
		matrix M(m_cols, m_rows);
		for(uint32_t i = 0; i < m_rows; i++)
		{
			for(uint32_t j = 0; j < m_cols; j++)
			{
				M.m_ptr[j + m_cols * i] = m_ref[i + m_rows * j];
			}
		}
		return M;
	}
	bool matrix::solve(matrix& x, const matrix& f) const
	{
		//data
		x = f;
		int32_t status;
		matrix M(*this);
		uint32_t* pivot = (uint32_t*) alloca(m_rows * sizeof(uint32_t));
		//check
		if(m_rows != m_cols)
		{
			fprintf(stderr, "Error: solve called on non-square matrix!\n");
			exit(EXIT_FAILURE);
		}
		if(m_cols != x.m_rows || x.m_rows != f.m_rows || x.m_cols != f.m_cols)
		{
			fprintf(stderr, "Error: solve called with incompatible matrices!\n");
			exit(EXIT_FAILURE);
		}
		//solve
		dgesv_(&m_rows, &x.m_cols, M.m_ptr, &m_rows, pivot, x.m_ptr, &m_rows, &status);
		//return
		return status == 0;
	}

	bool matrix::symmetric(double t) const
	{
		if(m_rows != m_cols)
		{
			return false;
		}
		else
		{
			double s = 0, v = 0;
			for(uint32_t i = 0; i < m_rows; i++)
			{
				s = fmax(s, m_ref[i + m_rows * i]);
				for(uint32_t j = i + 1; j < m_cols; j++)
				{
					s = fmax(s, m_ref[i + m_rows * j]);
					s = fmax(s, m_ref[j + m_rows * i]);
					v = fmax(v, fabs(m_ref[i + m_rows * j] - m_ref[j + m_rows * i]));
				}
			}
			return m_rows == 1 || v < t * s;
		}
	}

	mat3 matrix::span3(uint32_t r, uint32_t c) const
	{
		mat3 m;
		for(uint32_t i = 0; i < 3; i++)
		{
			for(uint32_t j = 0; j < 3; j++)
			{
				m[i + 3 * j] = m_ref[i + r + m_rows * (j + c)];
			}
		}
		return m;
	}
	math::span matrix::span(uint32_t r, uint32_t c, uint32_t n, uint32_t m)
	{
		return math::span(*this, r, c, n, m);
	}

	//svd
	bool matrix::svd(matrix& U, matrix& V, vector& s) const
	{
		//check
		if(s.m_rows != std::min(m_rows, m_cols))
		{
			fprintf(stderr, "\tError: SVD third parameter has incompatible dimensions!\n");
		}
		if(U.m_rows != U.m_cols || U.m_rows != m_rows)
		{
			fprintf(stderr, "\tError: SVD first parameter has incompatible dimensions!\n");
		}
		if(V.m_rows != V.m_cols || V.m_cols != m_cols)
		{
			fprintf(stderr, "\tError: SVD second parameter has incompatible dimensions!\n");
		}
		//data
		int32_t info;
		double query;
		matrix A(*this);
		uint32_t lwork = -1;
		const char* jobz = "A";
		int32_t* iwork = (int32_t*) alloca(8 * std::min(m_rows, m_cols) * sizeof(int32_t));
		//query
		dgesdd_(jobz, &m_rows, &m_cols, A.m_ptr, &m_rows, s.m_ptr, U.m_ptr, &m_rows, V.m_ptr, &m_cols, &query, &lwork, iwork, &info);
		//decompose
		lwork = (uint32_t) query;
		double* work = (double*) alloca(lwork * sizeof(double));
		dgesdd_(jobz, &m_rows, &m_cols, A.m_ptr, &m_rows, s.m_ptr, U.m_ptr, &m_rows, V.m_ptr, &m_cols, work, &lwork, iwork, &info);
		//return
		return info == 0;
	}

	//eigen
	bool matrix::eigen_sym(vector& v, matrix& P) const
	{
		//data
		double query, *work;
		int32_t status, lwork = -1;
		//query
		dsyev_("V", "L", &m_rows, nullptr, &m_rows, nullptr, &query, &lwork, &status);
		//eigen
		P = *this;
		lwork = (int32_t) query;
		work = new double[lwork];
		dsyev_("V", "L", &m_rows, P.m_ptr, &m_rows, v.m_ptr, work, &lwork, &status);
		//delete
		delete[] work;
		//return
		return status == 0;
	}
	bool matrix::eigen_sym(vector& v, matrix& P, const matrix& B) const
	{
		//data
		double query, *work;
		const uint32_t type = 1;
		int32_t status, lwork = -1;
		//query
		dsygv_(&type, "V", "L", &m_rows, nullptr, &m_rows, nullptr, &m_rows, nullptr, &query, &lwork, &status);
		//eigen
		P = *this;
		matrix Q = B;
		lwork = (int32_t) query;
		work = new double[lwork];
		dsygv_(&type, "V", "L", &m_rows, P.m_ptr, &m_rows, Q.m_ptr, &m_rows, v.m_ptr, work, &lwork, &status);
		//delete
		delete[] work;
		//return
		return status == 0;
	}

	//stats
	double matrix::mean(void) const
	{
		double m = 0;
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			m += m_ref[i];
		}
		return m / m_rows / m_cols;
	}
	double matrix::variance(void) const
	{
		double s = 0;
		const double m = mean();
		for(uint32_t i = 0; i < m_rows * m_cols; i++)
		{
			s += m_ref[i] * m_ref[i];
		}
		return s / m_rows / m_cols - m * m;
	}

	//static
	matrix matrix::eye(uint32_t rows, uint32_t cols)
	{
		return matrix(rows, cols).eye();
	}
	matrix matrix::zeros(uint32_t rows, uint32_t cols)
	{
		return matrix(rows, cols).zeros();
	}

	//friends
	matrix operator*(double s, const matrix& m)
	{
		return matrix(m) *= s;
	}
}