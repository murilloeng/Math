//std
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdint>

//Math
#include "Math/Test/inc/eigen.hpp"
#include "Math/inc/Linear/Eigen.hpp"
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Miscellaneous/util.hpp"

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
	//base
	for(uint32_t i = 0; i < order; i++)
	{
		for(uint32_t j = 0; j < order; j++)
		{
			A[i + order * j] = math::randu();
		}
	}
	//matrix
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

void tests::eigen::symmetric_std_full(void)
{
	//data
	math::Eigen eigen;
	const uint32_t order_max = 100;
	double A[order_max * order_max];
	//test
	eigen.data(0, A);
	srand(time(nullptr));
	for(uint32_t i = 1; i <= order_max; i++)
	{
		eigen.order(i);
		setup_symmetric_matrix(A, i);
		bool test = eigen.compute(true);
		for(uint32_t j = 0; j < i; j++)
		{
			const double w = eigen.eigenvalue(0, j);
			const double* z = eigen.eigenvector(0, j);
			test = test && fabs(math::Matrix(A, i, i).bilinear(z) - w) < 1e-5;
		}
		if(!test) break;
		printf("Test %3d: ok!\n", i);
	}
}
void tests::eigen::symmetric_gen_full(void)
{
	//data
	math::Eigen eigen;
	const uint32_t order_max = 100;
	double A[order_max * order_max];
	double B[order_max * order_max];
	//test
	eigen.data(0, A);
	eigen.data(1, B);
	srand(time(nullptr));
	for(uint32_t i = 1; i <= order_max; i++)
	{
		eigen.order(i);
		setup_symmetric_matrix(A, i);
		setup_symmetric_pd_matrix(B, i);
		bool test = eigen.compute(true);
		for(uint32_t j = 0; j < i; j++)
		{
			const double w = eigen.eigenvalue(0, j);
			const double* z = eigen.eigenvector(0, j);
			test = test && fabs(math::Matrix(A, i, i).bilinear(z) - w) < 1e-5;
			test = test && fabs(math::Matrix(B, i, i).bilinear(z) - 1) < 1e-5;
		}
		if(!test) break;
		printf("Test %3d: ok!\n", i);
	}
}
void tests::eigen::symmetric_std_partial(void)
{
	//data
	math::Eigen eigen;
	const uint32_t modes = 5;
	const uint32_t order_max = 100;
	double A[order_max * order_max];
	//test
	eigen.data(0, A);
	eigen.index_min(1);
	srand(time(nullptr));
	eigen.index_max(modes);
	eigen.type(math::Eigen::Type::Index);
	for(uint32_t i = modes; i <= order_max; i++)
	{
		eigen.order(i);
		setup_symmetric_matrix(A, i);
		bool test = eigen.compute(true);
		for(uint32_t j = 0; j < modes; j++)
		{
			const double w = eigen.eigenvalue(0, j);
			const double* z = eigen.eigenvector(0, j);
			test = test && fabs(math::Matrix(A, i, i).bilinear(z) - w) < 1e-5;
		}
		if(!test) break;
		printf("Test %3d: ok!\n", i);
	}
}
void tests::eigen::symmetric_gen_partial(void)
{
	//data
	math::Eigen eigen;
	const uint32_t modes = 5;
	const uint32_t order_max = 100;
	double A[order_max * order_max];
	double B[order_max * order_max];
	//test
	eigen.data(0, A);
	eigen.data(1, B);
	eigen.index_min(1);
	srand(time(nullptr));
	eigen.index_max(modes);
	eigen.type(math::Eigen::Type::Index);
	for(uint32_t i = modes; i <= order_max; i++)
	{
		eigen.order(i);
		setup_symmetric_matrix(A, i);
		setup_symmetric_pd_matrix(B, i);
		bool test = eigen.compute(true);
		for(uint32_t j = 0; j < modes; j++)
		{
			const double w = eigen.eigenvalue(0, j);
			const double* z = eigen.eigenvector(0, j);
			test = test && fabs(math::Matrix(A, i, i).bilinear(z) - w) < 1e-5;
			test = test && fabs(math::Matrix(B, i, i).bilinear(z) - 1) < 1e-5;
		}
		if(!test) break;
		printf("Test %3d: ok!\n", i);
	}
}