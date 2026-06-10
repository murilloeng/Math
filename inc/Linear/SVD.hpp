#pragma once

//std
#include <cstdint>

namespace math
{
	class SVD
	{
	public:
		//constructor
		SVD(void);

		//destructor
		~SVD(void);

		//data
		bool modes(bool);
		bool modes(void) const;

		uint32_t rows(uint32_t);
		uint32_t rows(void) const;

		uint32_t cols(uint32_t);
		uint32_t cols(void) const;

		const double* data(void) const;
		const double* data(const double*);

		//compute
		bool compute(void);

		//singular values
		double singular_value(uint32_t) const;
		const double* singular_values(void) const;

		//singular modes
		const double* singular_modes(uint32_t) const;
		const double* singular_modes(uint32_t, uint32_t) const;

	private:
		//setup
		void cleanup(void);
		void allocate(void);

		//data
		bool m_modes;
		uint32_t m_rows;
		uint32_t m_cols;
		const double* m_data;
		double* m_singular_values;
		double* m_singular_modes[2];

	};
}