#pragma once

//Math
#include "Math/inc/Eigen/SparseNon.hpp"

namespace math
{
	namespace eigen
	{
		class SparseNonGen : public SparseNon
		{
		public:
			//constructor
			SparseNonGen(void);
			SparseNonGen(uint32_t, uint32_t, uint32_t, const int32_t*, const int32_t*, const double*, const double*, double*, double*, double*);

			//destructor
			~SparseNonGen(void);

			//compute
			bool compute(void) override;
			void setup(uint32_t, uint32_t, uint32_t, const int32_t*, const int32_t*, const double*, const double*, double*, double*, double*);

		protected:
			//operation
			void operation(double*, double*) const override;

			//data
			double* m_U;
			double* m_sr;
			double* m_si;
			const double* m_A;
			const double* m_B;
			const math::Sparse* m_S;
		};
	}
}