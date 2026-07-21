#pragma once

//Math
#include "Math/inc/Eigen/SparseNon.hpp"

namespace math
{
	namespace eigen
	{
		class SparseNonStd : public SparseNon
		{
		public:
			//constructor
			SparseNonStd(void);
			SparseNonStd(uint32_t, uint32_t, uint32_t, const int32_t*, const int32_t*, const double*, double*, double*, double*);

			//destructor
			~SparseNonStd(void);

			//compute
			bool compute(void) override;
			void setup(uint32_t, uint32_t, uint32_t, const int32_t*, const int32_t*, const double*, double*, double*, double*);

		protected:
			//operation
			void operation(double*, double*) const override;

			//data
			double* m_U;
			double* m_sr;
			double* m_si;
			const double* m_A;
		};
	}
}