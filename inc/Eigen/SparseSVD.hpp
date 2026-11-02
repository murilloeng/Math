#pragma once

//Math
#include "Math/inc/Eigen/SparseSym.hpp"

namespace math
{
	namespace eigen
	{
		class SparseSVD : public SparseSym
		{
		public:
			//constructor
			SparseSVD(void);
			SparseSVD(uint32_t, uint32_t, uint32_t, uint32_t, const int32_t*, const int32_t*, const double*, double*, double*, double*);

			//destructor
			~SparseSVD(void) = default;

			//data
			uint32_t rows(void) const;
			uint32_t cols(void) const;

			//compute
			bool compute(void) override;
			void setup(uint32_t, uint32_t, uint32_t, uint32_t, const int32_t*, const int32_t*, const double*, double*, double*, double*);

		protected:
			//operation
			void operation(double*, double*) const override;

			//data
			uint32_t m_rows;
			double* m_s;
			double* m_U;
			double* m_V;
			double* m_w;
			const double* m_A;
		};
	}
}
