#pragma once

//std
#include <cstdint>

namespace math
{
	class SVD
	{
	public:
		//constructor
		SVD(const double*, uint32_t, uint32_t);
		SVD(const double*, uint32_t, uint32_t, double*, double*, double*);

		//destructor
		~SVD(void);

		//data
		bool modes(bool);
		bool modes(void) const;

		uint32_t rows(void) const;
		uint32_t cols(void) const;
		const double* data(void) const;
		const double* singular_values(void) const;
		const double* singular_vectors(uint32_t) const;

		//compute
		bool compute(void);

	private:
		//setup
		void cleanup(void);
		void allocate(void);

		//data
		bool m_own;
		bool m_modes;
		uint32_t m_rows;
		uint32_t m_cols;
		const double* m_A;
		double *m_s, *m_U, *m_V;

	};
}