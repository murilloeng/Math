//std
#include <cmath>
#include <ctime>
#include <cstring>
#include <cstdint>
#include <stdexcept>

//Suitesparse
#include "suitesparse/umfpack.h"

//Math
#include "Math/Test/inc/eigen.hpp"
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Linear/Sparse.hpp"
#include "Math/inc/Miscellaneous/util.hpp"

#include "Math/inc/Eigen/SVD.hpp"
#include "Math/inc/Eigen/DenseSymStd.hpp"
#include "Math/inc/Eigen/DenseSymGen.hpp"
#include "Math/inc/Eigen/DenseNonStd.hpp"
#include "Math/inc/Eigen/DenseNonGen.hpp"

#include "Math/inc/Eigen/SparseSymStd.hpp"
#include "Math/inc/Eigen/SparseSymGen.hpp"
#include "Math/inc/Eigen/SparseNonStd.hpp"
#include "Math/inc/Eigen/SparseNonGen.hpp"

static void random_dense_non(double* A, uint32_t order)
{
	for(uint32_t i = 0; i < order; i++)
	{
		for(uint32_t j = 0; j < order; j++)
		{
			A[i + order * j] = math::randu();
		}
	}
}
static void random_dense_sym(double* A, uint32_t order)
{
	for(uint32_t i = 0; i < order; i++)
	{
		for(uint32_t j = i; j < order; j++)
		{
			A[i + order * j] = A[j + order * i] = math::randu();
		}
	}
}
static void random_dense_sym_pd(double* B, uint32_t order)
{
	//data
	double* A = new double[order * order];
	//matrix
	random_dense_non(A, order);
	for(uint32_t i = 0; i < order; i++)
	{
		for(uint32_t j = 0; j < order; j++)
		{
			B[i + order * j] = 0;
			for(uint32_t k = 0; k < order; k++)
			{
				B[i + order * j] += A[k + order * i] * A[k + order * j];
			}
		}
	}
	//return
	delete[] A;
}

static void random_sparse_sym(double*& Am, int32_t*& rm, int32_t*& cm, uint32_t n)
{
	//data
	const double v_min = -1.00e+00;
	const double v_max = +1.00e+00;
	const double density = 5.00e-02;
	//triplets
	srand(time(nullptr));
	std::vector<double> At;
	std::vector<int32_t> rt, ct;
	for(uint32_t i = 0; i < n; i++)
	{
		for(uint32_t j = i; j < n; j++)
		{
			if(i == j || math::randu() < density)
			{
				//data
				const double v = math::randu(v_min, v_max);
				//append
				rt.push_back(i);
				ct.push_back(j);
				At.push_back(v);
				if(i != j) rt.push_back(j);
				if(i != j) ct.push_back(i);
				if(i != j) At.push_back(v);
			}
		}
	}
	//sparse
	cm = new int32_t[n + 1];
	rm = new int32_t[At.size()];
	Am = new double[At.size()];
	umfpack_di_triplet_to_col(n, n, At.size(), rt.data(), ct.data(), At.data(), cm, rm, Am, nullptr);
}
static void random_sparse_non(double*& Am, int32_t*& rm, int32_t*& cm, uint32_t n)
{
	//data
	const double v_min = -1.00e+00;
	const double v_max = +1.00e+00;
	const double density = 5.00e-02;
	//triplets
	srand(time(nullptr));
	std::vector<double> At;
	std::vector<int32_t> rt, ct;
	for(uint32_t i = 0; i < n; i++)
	{
		for(uint32_t j = i; j < n; j++)
		{
			if(math::randu() < density)
			{
				rt.push_back(i);
				ct.push_back(j);
				At.push_back(math::randu(v_min, v_max));
			}
		}
	}
	//sparse
	cm = new int32_t[n + 1];
	rm = new int32_t[At.size()];
	Am = new double[At.size()];
	umfpack_di_triplet_to_col(n, n, At.size(), rt.data(), ct.data(), At.data(), cm, rm, Am, nullptr);
}
static void random_sparse_sym_pd(double*& Bm, const int32_t* rm, const int32_t* cm, uint32_t n)
{
	//data
	const double v_min = -1.00e+00;
	const double v_max = +1.00e+00;
	//values
	srand(time(nullptr));
	Bm = new double[cm[n]];
	for(int32_t i = 0; i < cm[n]; i++)
	{
		Bm[i] = math::randu(v_min, v_max);
	}
	//symmetry
	math::Sparse Bs(Bm, rm, cm, n, n);
	for(uint32_t i = 0; i < n; i++)
	{
		for(int32_t j = cm[i]; j < cm[i + 1]; j++)
		{
			if(rm[j] < int32_t(i))
			{
				int32_t k;
				for(k = cm[rm[j]]; k < cm[rm[j] + 1]; k++) if(rm[k] == int32_t(i)) break;
				const double v1 = Bm[j];
				const double v2 = Bm[k];
				Bm[j] = Bm[k] = (v1 + v2) / 2;
			}
		}
	}
	//diagonal
	for(uint32_t i = 0; i < n; i++)
	{
		double d = 0;
		for(int32_t j = cm[i]; j < cm[i + 1]; j++)
		{
			if(rm[j] != int32_t(i)) d += fabs(Bm[j]);
		}
		for(int32_t j = cm[i]; j < cm[i + 1]; j++)
		{
			if(rm[j] == int32_t(i)) Bm[j] = d + 1;
		}
	}
}

void tests::eigen::dense_svd(void)
{
	//data
	srand(time(nullptr));
	const uint32_t order_max = 100;
	//test
	double s[order_max];
	double A[order_max * order_max];
	double B[order_max * order_max];
	double U[order_max * order_max];
	double V[order_max * order_max];
	for(uint32_t i = 1; i <= order_max; i++)
	{
		random_dense_non(A, i);
		memcpy(B, A, i * i * sizeof(double));
		bool test = math::eigen::SVD(A, i, i, s, U, V).compute();
		for(uint32_t j = 0; j < i; j++)
		{
			const double sj = s[j];
			const double* u = U + i * j;
			const double* v = V + i * j;
			test = test && fabs(math::Matrix(B, i, i).bilinear(u, v) - sj) < 1e-5;
		}
		if(!test) throw std::runtime_error("Error");
		printf("Test eigen dense svd %3d: ok!\n", i);
	}
}
void tests::eigen::dense_sym_std_full(void)
{
	//data
	srand(time(nullptr));
	const uint32_t order_max = 100;
	//test
	double s[order_max];
	double U[order_max * order_max];
	double A[order_max * order_max];
	double B[order_max * order_max];
	for(uint32_t i = 1; i <= order_max; i++)
	{
		random_dense_sym(A, i);
		memcpy(B, A, i * i * sizeof(double));
		bool test = math::eigen::DenseSymStd(i, A, s, U).compute();
		for(uint32_t j = 0; j < i; j++)
		{
			const double w = s[j];
			const double* z = U + i * j;
			test = test && fabs(math::Matrix(B, i, i).bilinear(z) - w) < 1e-5;
		}
		if(!test) throw std::runtime_error("Error");
		printf("Test eigen dense sym std full %3d: ok!\n", i);
	}
}
void tests::eigen::dense_sym_gen_full(void)
{
	//data
	srand(time(nullptr));
	const uint32_t order_max = 100;
	//test
	double s[order_max];
	double U[order_max * order_max];
	double A[order_max * order_max];
	double B[order_max * order_max];
	double C[order_max * order_max];
	double D[order_max * order_max];
	for(uint32_t i = 1; i <= order_max; i++)
	{
		random_dense_sym(A, i);
		random_dense_sym_pd(B, i);
		memcpy(C, A, i * i * sizeof(double));
		memcpy(D, B, i * i * sizeof(double));
		bool test = math::eigen::DenseSymGen(i, A, B, s, U).compute();
		for(uint32_t j = 0; j < i; j++)
		{
			const double w = s[j];
			const double* z = U + i * j;
			test = test && fabs(math::Matrix(C, i, i).bilinear(z) - w) < 1e-5;
			test = test && fabs(math::Matrix(D, i, i).bilinear(z) - 1) < 1e-5;
		}
		if(!test) throw std::runtime_error("Error");
		printf("Test eigen dense sym gen full %3d: ok!\n", i);
	}
}
void tests::eigen::dense_sym_std_part(void)
{
	//data
	srand(time(nullptr));
	const uint32_t modes = 5;
	const uint32_t order_max = 100;
	//test
	double s[modes];
	double U[modes * order_max];
	double A[order_max * order_max];
	double B[order_max * order_max];
	for(uint32_t i = modes; i <= order_max; i++)
	{
		random_dense_sym(A, i);
		memcpy(B, A, i * i * sizeof(double));
		bool test = math::eigen::DenseSymStd(i, A, s, U, 0, modes - 1).compute();
		for(uint32_t j = 0; j < modes; j++)
		{
			const double w = s[j];
			const double* z = U + i * j;
			test = test && fabs(math::Matrix(B, i, i).bilinear(z) - w) < 1e-5;
		}
		if(!test) throw std::runtime_error("Error");
		printf("Test eigen dense sym std part %3d: ok!\n", i);
	}
}
void tests::eigen::dense_sym_gen_part(void)
{
	//data
	srand(time(nullptr));
	const uint32_t modes = 5;
	const uint32_t order_max = 100;
	//test
	double s[modes];
	double U[modes * order_max];
	double A[order_max * order_max];
	double B[order_max * order_max];
	double C[order_max * order_max];
	double D[order_max * order_max];
	for(uint32_t i = modes; i <= order_max; i++)
	{
		random_dense_sym(A, i);
		random_dense_sym_pd(B, i);
		memcpy(C, A, i * i * sizeof(double));
		memcpy(D, B, i * i * sizeof(double));
		bool test = math::eigen::DenseSymGen(i, A, B, s, U, 0, modes - 1).compute();
		for(uint32_t j = 0; j < modes; j++)
		{
			const double w = s[j];
			const double* z = U + i * j;
			test = test && fabs(math::Matrix(C, i, i).bilinear(z) - w) < 1e-5;
			test = test && fabs(math::Matrix(D, i, i).bilinear(z) - 1) < 1e-5;
		}
		if(!test) throw std::runtime_error("Error");
		printf("Test eigen dense sym gen part %3d: ok!\n", i);
	}
}
void tests::eigen::dense_non_std_full(void)
{
	//data
	srand(time(nullptr));
	const uint32_t order_max = 100;
	//test
	double sr[order_max];
	double si[order_max];
	double U[order_max * order_max];
	double V[order_max * order_max];
	double A[order_max * order_max];
	double B[order_max * order_max];
	for(uint32_t i = 1; i <= order_max; i++)
	{
		math::Matrix Am(B, i, i);
		random_dense_non(A, i);
		memcpy(B, A, i * i * sizeof(double));
		bool test = math::eigen::DenseNonStd(i, A, sr, si, U, V).compute();
		for(uint32_t j = 0; j < i; j++)
		{
			const double wr = sr[j];
			const double wi = si[j];
			if(wi == 0)
			{
				const double* zr = U + i * j;
				test = test && fabs(Am.bilinear(zr, zr) - wr) < 1e-5;
			}
			else if(wi > 0)
			{
				math::Vector zr(U + i * (j + 0), i);
				math::Vector zi(U + i * (j + 1), i);
				test = test && fabs((Am.bilinear(zr, zr) + Am.bilinear(zi, zi)) / (zr.inner(zr) + zi.inner(zi)) - wr) < 1e-5;
				test = test && fabs((Am.bilinear(zr, zi) - Am.bilinear(zi, zr)) / (zr.inner(zr) + zi.inner(zi)) - wi) < 1e-5;
			}
			else if(wi < 0)
			{
				math::Vector zr(U + i * (j - 1), i);
				math::Vector zi(U + i * (j + 0), i);
				test = test && fabs((Am.bilinear(zr, zr) + Am.bilinear(zi, zi)) / (zr.inner(zr) + zi.inner(zi)) - wr) < 1e-5;
				test = test && fabs((Am.bilinear(zi, zr) - Am.bilinear(zr, zi)) / (zr.inner(zr) + zi.inner(zi)) - wi) < 1e-5;
			}
		}
		if(!test) throw std::runtime_error("Error");
		printf("Test eigen dense non std full %3d: ok!\n", i);
	}
}
void tests::eigen::dense_non_gen_full(void)
{
	//data
	srand(time(nullptr));
	const uint32_t order_max = 100;
	//test
	double sr[order_max];
	double si[order_max];
	double U[order_max * order_max];
	double V[order_max * order_max];
	double A[order_max * order_max];
	double B[order_max * order_max];
	double C[order_max * order_max];
	double D[order_max * order_max];
	for(uint32_t i = 1; i <= order_max; i++)
	{
		math::Matrix Am(C, i, i);
		math::Matrix Bm(D, i, i);
		random_dense_non(A, i);
		random_dense_non(B, i);
		memcpy(C, A, i * i * sizeof(double));
		memcpy(D, B, i * i * sizeof(double));
		bool test = math::eigen::DenseNonGen(i, A, B, sr, si, U, V).compute();
		for(uint32_t j = 0; j < i; j++)
		{
			const double wr = sr[j];
			const double wi = si[j];
			if(wi == 0)
			{
				const double* zr = U + i * j;
				test = test && fabs(Am.bilinear(zr, zr) - wr * Bm.bilinear(zr, zr)) < 1e-5;
			}
			else if(wi > 0)
			{
				math::Vector zr(U + i * (j + 0), i);
				math::Vector zi(U + i * (j + 1), i);
				const double Au = Am.bilinear(zr, zr) + Am.bilinear(zi, zi);
				const double Ac = Am.bilinear(zr, zi) - Am.bilinear(zi, zr);
				const double Bu = Bm.bilinear(zr, zr) + Bm.bilinear(zi, zi);
				const double Bc = Bm.bilinear(zr, zi) - Bm.bilinear(zi, zr);
				test = test && fabs(wr * Bu - wi * Bc - Au) < 1e-5 && fabs(wr * Bc + wi * Bu - Ac) < 1e-5;
			}
			else if(wi < 0)
			{
				math::Vector zr(U + i * (j - 1), i);
				math::Vector zi(U + i * (j + 0), i);
				const double Au = Am.bilinear(zr, zr) + Am.bilinear(zi, zi);
				const double Ac = Am.bilinear(zr, zi) - Am.bilinear(zi, zr);
				const double Bu = Bm.bilinear(zr, zr) + Bm.bilinear(zi, zi);
				const double Bc = Bm.bilinear(zr, zi) - Bm.bilinear(zi, zr);
				test = test && fabs(wr * Bu + wi * Bc - Au) < 1e-5 && fabs(wr * Bc - wi * Bu - Ac) < 1e-5;
			}
		}
		if(!test) throw std::runtime_error("Error");
		printf("Test eigen dense non gen full %3d: ok!\n", i);
	}
}

void tests::eigen::sparse_svd(void)
{
	return;
}
void tests::eigen::sparse_sym_std_part(void)
{
	int32_t *rm, *cm;
	const uint32_t n = 500;
	const uint32_t nev = 5;
	const uint32_t ncv = 15;
	double *Am, s[nev], U[n * nev];
	for(uint32_t i = 100; i <= n; i += 10)
	{
		random_sparse_sym(Am, rm, cm, i);
		bool test = math::eigen::SparseSymStd(i, nev, ncv, rm, cm, Am, s, U).compute();
		for(uint32_t j = 0; j < nev; j++)
		{
			test = test && fabs(math::Sparse(Am, rm, cm, i, i).bilinear(U + i * j) - s[j]) < 1e-5;
		}
		delete[] Am;
		delete[] rm;
		delete[] cm;
		if(!test) throw std::runtime_error("Error");
		printf("Test eigen sparse sym std part %3d: ok!\n", i);
	}
}
void tests::eigen::sparse_sym_gen_part(void)
{
	int32_t *rm, *cm;
	const uint32_t n = 500;
	const uint32_t nev = 5;
	const uint32_t ncv = 15;
	double *Am, *Bm, s[nev], U[n * nev];
	for(uint32_t i = 100; i <= n; i += 10)
	{
		random_sparse_sym(Am, rm, cm, i);
		random_sparse_sym_pd(Bm, rm, cm, i);
		bool test = math::eigen::SparseSymGen(i, nev, ncv, rm, cm, Am, Bm, s, U).compute();
		for(uint32_t j = 0; j < nev; j++)
		{
			test = test && fabs(math::Sparse(Bm, rm, cm, i, i).bilinear(U + i * j) - 1) < 1e-5;
			test = test && fabs(math::Sparse(Am, rm, cm, i, i).bilinear(U + i * j) - s[j]) < 1e-5;
		}
		delete[] Am;
		delete[] Bm;
		delete[] rm;
		delete[] cm;
		if(!test) throw std::runtime_error("Error");
		printf("Test eigen sparse sym gen part %3d: ok!\n", i);
	}
}
void tests::eigen::sparse_non_std_part(void)
{
	int32_t *rm, *cm;
	const uint32_t n = 500;
	const uint32_t nev = 5;
	const uint32_t ncv = 50;
	double *Am, sr[nev], si[nev];
	double* U = new double[n * nev];
	for(uint32_t i = 100; i <= n; i += 10)
	{
		random_sparse_non(Am, rm, cm, i);
		bool test = math::eigen::SparseNonStd(i, nev, ncv, rm, cm, Am, sr, si, U).compute();
		for(uint32_t j = 0; j < nev; j++)
		{
			const double wr = sr[j];
			const double wi = si[j];
			const math::Sparse As(Am, rm, cm, i, i);
			if(wi == 0)
			{
				const double* zr = U + i * j;
				test = test && fabs(As.bilinear(zr, zr) - wr) < 1e-5;
			}
			else if(wi > 0)
			{
				math::Vector zr(U + i * (j + 0), i);
				math::Vector zi(U + i * (j + 1), i);
				test = test && fabs((As.bilinear(zr, zr) + As.bilinear(zi, zi)) / (zr.inner(zr) + zi.inner(zi)) - wr) < 1e-5;
				test = test && fabs((As.bilinear(zr, zi) - As.bilinear(zi, zr)) / (zr.inner(zr) + zi.inner(zi)) - wi) < 1e-5;
			}
			else if(wi < 0)
			{
				math::Vector zr(U + i * (j - 1), i);
				math::Vector zi(U + i * (j + 0), i);
				test = test && fabs((As.bilinear(zr, zr) + As.bilinear(zi, zi)) / (zr.inner(zr) + zi.inner(zi)) - wr) < 1e-5;
				test = test && fabs((As.bilinear(zi, zr) - As.bilinear(zr, zi)) / (zr.inner(zr) + zi.inner(zi)) - wi) < 1e-5;
			}
		}
		delete[] Am;
		delete[] rm;
		delete[] cm;
		if(!test) throw std::runtime_error("Error");
		printf("Test eigen sparse non std part %3d: ok!\n", i);
	}
	delete[] U;
}
void tests::eigen::sparse_non_gen_part(void)
{
	int32_t *rm, *cm;
	const uint32_t n = 500;
	const uint32_t nev = 5;
	const uint32_t ncv = 50;
	double* U = new double[n * nev];
	double *Am, *Bm, sr[nev], si[nev];
	for(uint32_t i = 100; i <= n; i += 10)
	{
		random_sparse_non(Am, rm, cm, i);
		random_sparse_sym_pd(Bm, rm, cm, i);
		bool test = math::eigen::SparseNonGen(i, nev, ncv, rm, cm, Am, Bm, sr, si, U).compute();
		for(uint32_t j = 0; j < nev; j++)
		{
			const double wr = sr[j];
			const double wi = si[j];
			const math::Sparse As(Am, rm, cm, i, i);
			const math::Sparse Bs(Bm, rm, cm, i, i);
			if(wi == 0)
			{
				const double* zr = U + i * j;
				test = test && fabs(As.bilinear(zr, zr) - wr * Bs.bilinear(zr, zr)) < 1e-5;
			}
			else if(wi > 0)
			{
				math::Vector zr(U + i * (j + 0), i);
				math::Vector zi(U + i * (j + 1), i);
				const double Au = As.bilinear(zr, zr) + As.bilinear(zi, zi);
				const double Ac = As.bilinear(zr, zi) - As.bilinear(zi, zr);
				const double Bu = Bs.bilinear(zr, zr) + Bs.bilinear(zi, zi);
				const double Bc = Bs.bilinear(zr, zi) - Bs.bilinear(zi, zr);
				test = test && fabs(wr * Bu - wi * Bc - Au) < 1e-5 && fabs(wr * Bc + wi * Bu - Ac) < 1e-5;
			}
			else if(wi < 0)
			{
				math::Vector zr(U + i * (j - 1), i);
				math::Vector zi(U + i * (j + 0), i);
				const double Au = As.bilinear(zr, zr) + As.bilinear(zi, zi);
				const double Ac = As.bilinear(zr, zi) - As.bilinear(zi, zr);
				const double Bu = Bs.bilinear(zr, zr) + Bs.bilinear(zi, zi);
				const double Bc = Bs.bilinear(zr, zi) - Bs.bilinear(zi, zr);
				test = test && fabs(wr * Bu + wi * Bc - Au) < 1e-5 && fabs(wr * Bc - wi * Bu - Ac) < 1e-5;
			}
		}
		delete[] Am;
		delete[] rm;
		delete[] cm;
		if(!test) throw std::runtime_error("Error");
		printf("Test eigen sparse non gen part %3d: ok!\n", i);
	}
	delete[] U;
}