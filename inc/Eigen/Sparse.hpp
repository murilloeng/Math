#pragma once

//Math
#include "Math/inc/Eigen/Eigen.hpp"

namespace math
{
	namespace eigen
	{
		class Sparse : public Eigen
		{
		public:
			//constructor
			Sparse(uint32_t, uint32_t, uint32_t, const int32_t*, const int32_t*);

			//destructor
			virtual ~Sparse(void);

			//data
			uint32_t mode(uint32_t);
			uint32_t mode(void) const;
			
			uint32_t modes(uint32_t);
			uint32_t modes(void) const;
			
			uint32_t vectors(uint32_t);
			uint32_t vectors(void) const;
			
			uint32_t iteration_max(uint32_t);
			uint32_t iteration_max(void) const;

			double shift(double);
			double shift(void) const;

			const char* type(void) const;
			const char* type(const char*);

			const int32_t* rows_map(void) const;
			const int32_t* rows_map(const int32_t*);

			const int32_t* cols_map(void) const;
			const int32_t* cols_map(const int32_t*);

		protected:
			//operation
			virtual void operation(double*, double*) const = 0;

			//data
			int32_t m_ido;
			uint32_t m_mode;
			uint32_t m_modes;
			uint32_t m_vectors;
			uint32_t m_iteration_max;

			char m_type[3];
			double m_shift;
			const int32_t* m_rows_map;
			const int32_t* m_cols_map;
		};
	}
}