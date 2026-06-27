//std
#include <cmath>
#include <ctime>
#include <cstring>
#include <cstdint>
#include <stdexcept>

//Math
#include "Math/Test/inc/eigen.hpp"
#include "Math/inc/Linear/SVD.hpp"
#include "Math/inc/Linear/Eigen.hpp"
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Miscellaneous/util.hpp"

static void setup_random_matrix(double* A, uint32_t order)
{
	for(uint32_t i = 0; i < order; i++)
	{
		for(uint32_t j = 0; j < order; j++)
		{
			A[i + order * j] = math::randu();
		}
	}
}
static void setup_symmetric_matrix(double* A, uint32_t order)
{
	for(uint32_t i = 0; i < order; i++)
	{
		for(uint32_t j = i; j < order; j++)
		{
			A[i + order * j] = A[j + order * i] = math::randu();
		}
	}
}
static void setup_symmetric_pd_matrix(double* B, uint32_t order)
{
	//data
	double* A = new double[order * order];
	//matrix
	setup_random_matrix(A, order);
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

void tests::eigen::non_symmetric_std(void)
{
	//data
	math::Eigen eigen;
	srand(time(nullptr));
	const uint32_t order_max = 100;
	double A[order_max * order_max];
	//test
	eigen.data(0, A);
	eigen.symmetric(false);
	for(uint32_t i = 1; i <= order_max; i++)
	{
		eigen.order(i);
		math::Matrix Am(A, i, i);
		setup_random_matrix(A, i);
		bool test = eigen.compute();
		for(uint32_t j = 0; j < i; j++)
		{
			const double wr = eigen.eigenvalue(0, j);
			const double wi = eigen.eigenvalue(1, j);
			if(wi == 0)
			{
				const double* zr = eigen.eigenvector(0, j);
				test = test && fabs(Am.bilinear(zr, zr) - wr) < 1e-5;
			}
			else if(wi > 0)
			{
				math::Vector zr(eigen.eigenvector(0, j + 0), i);
				math::Vector zi(eigen.eigenvector(0, j + 1), i);
				test = test && fabs((Am.bilinear(zr, zr) + Am.bilinear(zi, zi)) / (zr.inner(zr) + zi.inner(zi)) - wr) < 1e-5;
				test = test && fabs((Am.bilinear(zr, zi) - Am.bilinear(zi, zr)) / (zr.inner(zr) + zi.inner(zi)) - wi) < 1e-5;
			}
			else if(wi < 0)
			{
				math::Vector zr(eigen.eigenvector(0, j - 1), i);
				math::Vector zi(eigen.eigenvector(0, j + 0), i);
				test = test && fabs((Am.bilinear(zr, zr) + Am.bilinear(zi, zi)) / (zr.inner(zr) + zi.inner(zi)) - wr) < 1e-5;
				test = test && fabs((Am.bilinear(zi, zr) - Am.bilinear(zr, zi)) / (zr.inner(zr) + zi.inner(zi)) - wi) < 1e-5;
			}
		}
		if(!test) throw std::runtime_error("Error");
		printf("Test non symmetric std %3d: ok!\n", i);
	}
}
void tests::eigen::non_symmetric_gen(void)
{
	//data
	math::Eigen eigen;
	srand(time(nullptr));
	const uint32_t order_max = 100;
	double A[order_max * order_max];
	double B[order_max * order_max];
	//test
	eigen.data(0, A);
	eigen.data(1, B);
	eigen.symmetric(false);
	for(uint32_t i = 1; i <= order_max; i++)
	{
		eigen.order(i);
		math::Matrix Am(A, i, i);
		math::Matrix Bm(B, i, i);
		setup_random_matrix(A, i);
		setup_random_matrix(B, i);
		bool test = eigen.compute();
		for(uint32_t j = 0; j < i; j++)
		{
			const double wr = eigen.eigenvalue(0, j);
			const double wi = eigen.eigenvalue(1, j);
			if(wi == 0)
			{
				const double* zr = eigen.eigenvector(0, j);
				test = test && fabs(Am.bilinear(zr, zr) - wr * Bm.bilinear(zr, zr)) < 1e-5;
			}
			else if(wi > 0)
			{
				math::Vector zr(eigen.eigenvector(0, j + 0), i);
				math::Vector zi(eigen.eigenvector(0, j + 1), i);
				const double Au = Am.bilinear(zr, zr) + Am.bilinear(zi, zi);
				const double Ac = Am.bilinear(zr, zi) - Am.bilinear(zi, zr);
				const double Bu = Bm.bilinear(zr, zr) + Bm.bilinear(zi, zi);
				const double Bc = Bm.bilinear(zr, zi) - Bm.bilinear(zi, zr);
				test = test && fabs(wr * Bu - wi * Bc - Au) < 1e-5 && fabs(wr * Bc + wi * Bu - Ac) < 1e-5;
			}
			else if(wi < 0)
			{
				math::Vector zr(eigen.eigenvector(0, j - 1), i);
				math::Vector zi(eigen.eigenvector(0, j + 0), i);
				const double Au = Am.bilinear(zr, zr) + Am.bilinear(zi, zi);
				const double Ac = Am.bilinear(zr, zi) - Am.bilinear(zi, zr);
				const double Bu = Bm.bilinear(zr, zr) + Bm.bilinear(zi, zi);
				const double Bc = Bm.bilinear(zr, zi) - Bm.bilinear(zi, zr);
				test = test && fabs(wr * Bu + wi * Bc - Au) < 1e-5 && fabs(wr * Bc - wi * Bu - Ac) < 1e-5;
			}
		}
		if(!test) throw std::runtime_error("Error");
		printf("Test non symmetric gen %3d: ok!\n", i);
	}
}
void tests::eigen::symmetric_std_full(void)
{
	//data
	math::Eigen eigen;
	srand(time(nullptr));
	const uint32_t order_max = 100;
	double A[order_max * order_max];
	//test
	eigen.data(0, A);
	for(uint32_t i = 1; i <= order_max; i++)
	{
		eigen.order(i);
		setup_symmetric_matrix(A, i);
		bool test = eigen.compute();
		for(uint32_t j = 0; j < i; j++)
		{
			const double w = eigen.eigenvalue(0, j);
			const double* z = eigen.eigenvector(0, j);
			test = test && fabs(math::Matrix(A, i, i).bilinear(z) - w) < 1e-5;
		}
		if(!test) throw std::runtime_error("Error");
		printf("Test symmetric std full %3d: ok!\n", i);
	}
}
void tests::eigen::symmetric_gen_full(void)
{
	//data
	math::Eigen eigen;
	srand(time(nullptr));
	const uint32_t order_max = 100;
	double A[order_max * order_max];
	double B[order_max * order_max];
	//test
	eigen.data(0, A);
	eigen.data(1, B);
	for(uint32_t i = 1; i <= order_max; i++)
	{
		eigen.order(i);
		setup_symmetric_matrix(A, i);
		setup_symmetric_pd_matrix(B, i);
		bool test = eigen.compute();
		for(uint32_t j = 0; j < i; j++)
		{
			const double w = eigen.eigenvalue(0, j);
			const double* z = eigen.eigenvector(0, j);
			test = test && fabs(math::Matrix(A, i, i).bilinear(z) - w) < 1e-5;
			test = test && fabs(math::Matrix(B, i, i).bilinear(z) - 1) < 1e-5;
		}
		if(!test) throw std::runtime_error("Error");
		printf("Test symmetric gen full %3d: ok!\n", i);
	}
}
void tests::eigen::symmetric_std_partial(void)
{
	//data
	math::Eigen eigen;
	srand(time(nullptr));
	const uint32_t modes = 5;
	const uint32_t order_max = 100;
	double A[order_max * order_max];
	//test
	eigen.data(0, A);
	eigen.index_min(0);
	eigen.index_max(modes - 1);
	eigen.type(math::Eigen::Type::Index);
	for(uint32_t i = modes; i <= order_max; i++)
	{
		eigen.order(i);
		setup_symmetric_matrix(A, i);
		bool test = eigen.compute();
		for(uint32_t j = 0; j < modes; j++)
		{
			const double w = eigen.eigenvalue(0, j);
			const double* z = eigen.eigenvector(0, j);
			test = test && fabs(math::Matrix(A, i, i).bilinear(z) - w) < 1e-5;
		}
		if(!test) throw std::runtime_error("Error");
		printf("Test symmetric std partial %3d: ok!\n", i);
	}
}
void tests::eigen::symmetric_gen_partial(void)
{
	//data
	math::Eigen eigen;
	srand(time(nullptr));
	const uint32_t modes = 5;
	const uint32_t order_max = 100;
	double A[order_max * order_max];
	double B[order_max * order_max];
	//test
	eigen.data(0, A);
	eigen.data(1, B);
	eigen.index_min(0);
	eigen.index_max(modes - 1);
	eigen.type(math::Eigen::Type::Index);
	for(uint32_t i = modes; i <= order_max; i++)
	{
		eigen.order(i);
		setup_symmetric_matrix(A, i);
		setup_symmetric_pd_matrix(B, i);
		bool test = eigen.compute();
		for(uint32_t j = 0; j < modes; j++)
		{
			const double w = eigen.eigenvalue(0, j);
			const double* z = eigen.eigenvector(0, j);
			test = test && fabs(math::Matrix(A, i, i).bilinear(z) - w) < 1e-5;
			test = test && fabs(math::Matrix(B, i, i).bilinear(z) - 1) < 1e-5;
		}
		if(!test) throw std::runtime_error("Error");
		printf("Test symmetric gen partial %3d: ok!\n", i);
	}
}
void tests::eigen::singular_value_decomposition(void)
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
		setup_random_matrix(A, i);
		memcpy(B, A, i * i * sizeof(double));
		bool test = math::SVD(A, i, i, s, U, V).compute();
		for(uint32_t j = 0; j < i; j++)
		{
			const double sj = s[j];
			const double* u = U + i * j;
			const double* v = V + i * j;
			test = test && fabs(math::Matrix(B, i, i).bilinear(u, v) - sj) < 1e-5;
		}
		if(!test) throw std::runtime_error("Error");
		printf("Test singular value decomposition %3d: ok!\n", i);
	}
}