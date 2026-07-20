#pragma once

//Math
#include "Math/inc/Linear/Sparse.hpp"
#include "Math/inc/Eigen/SparseSym.hpp"

namespace math
{
	namespace eigen
	{
		class SparseSymGen : public SparseSym
		{
		public:
			//constructor
			SparseSymGen(void);
			SparseSymGen(uint32_t, uint32_t, uint32_t, const int32_t*, const int32_t*, const double*, const double*, double*, double*);

			//destructor
			~SparseSymGen(void);

			//compute
			bool compute(void) override;
			void setup(uint32_t, uint32_t, uint32_t, const int32_t*, const int32_t*, const double*, const double*, double*, double*);

		protected:
			//operation
			void operation(double*, double*) const override;

			//data
			double* m_s;
			double* m_U;
			const double* m_A;
			const double* m_B;
			const math::Sparse* m_S;
		};
	}
}