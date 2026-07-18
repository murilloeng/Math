#pragma once

//Math
#include "Math/inc/Eigen/EigenDenseStd.hpp"

namespace math
{
	namespace eigen
	{
		class EigenDenseStdSym : public EigenDenseStd
		{
		public:
			//constructor
			EigenDenseStdSym(uint32_t, double*, double*, double*);

			//destructor
			virtual ~EigenDenseStdSym(void);

			//data
			char range(char);
			char range(void) const;

			uint32_t modes(void) const;

			double value_min(double);
			double value_min(void) const;

			double value_max(double);
			double value_max(void) const;

			uint32_t index_min(uint32_t);
			uint32_t index_min(void) const;

			uint32_t index_max(uint32_t);
			uint32_t index_max(void) const;

			//compute
			bool compute(void) override;

		protected:
			//compute
			bool compute_full(void);
			bool compute_partial(void);

			//data
			double* m_s;
			double* m_U;
			char m_range;
			uint32_t m_modes;
			double m_value_min;
			double m_value_max;
			uint32_t m_index_min;
			uint32_t m_index_max;
		};
	}
}