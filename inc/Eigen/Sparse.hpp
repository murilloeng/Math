#pragma once

//Math
#include "Math/inc/Eigen/Base.hpp"

namespace math
{
	namespace eigen
	{
		class Sparse : public Base
		{
		public:
			//constructor
			Sparse(uint32_t, uint32_t, uint32_t, const int32_t*, const int32_t*);

			//destructor
			virtual ~Sparse(void);

			//data
			uint32_t modes(uint32_t);
			uint32_t modes(void) const;

			uint32_t vectors(uint32_t);
			uint32_t vectors(void) const;

			const int32_t* rows_map(void) const;
			const int32_t* rows_map(const int32_t*);

			const int32_t* cols_map(void) const;
			const int32_t* cols_map(const int32_t*);

		protected:
			//data
			uint32_t m_modes;
			uint32_t m_vectors;

			const int32_t* m_rows_map;
			const int32_t* m_cols_map;
		};
	}
}