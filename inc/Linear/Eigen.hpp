#pragma once

//std
#include <cstdint>

namespace math
{
	class Eigen
	{
	public:
		//constructor
		Eigen(double*, uint32_t, double*, double*);
		Eigen(double*, uint32_t, double*, double*, double*, double*);
		Eigen(double*, uint32_t, double*, double*, uint32_t, uint32_t);

		Eigen(double*, double*, uint32_t, double*, double*);
		Eigen(double*, double*, uint32_t, double*, double*, double*, double*);
		Eigen(double*, double*, uint32_t, double*, double*, uint32_t, uint32_t);

		//destructor
		~Eigen(void);

		//data
		bool full(void) const;

		bool modes(void) const;

		bool symmetric(void) const;

		uint32_t order(void) const;

		double* data(uint32_t) const;

		uint32_t index_min(void) const;
		uint32_t index_max(void) const;

		//compute
		bool compute(void);

		//eigenvalues
		const double* eigenvalues(uint32_t) const;
		double eigenvalue(uint32_t, uint32_t) const;

		//eigenvectors
		const double* eigenvectors(uint32_t) const;
		const double* eigenvector(uint32_t, uint32_t) const;

	private:
		//compute
		bool compute_non_symmetric_std(void);
		bool compute_non_symmetric_gen(void);
		bool compute_symmetric_std_full(void);
		bool compute_symmetric_gen_full(void);
		bool compute_symmetric_std_partial(void);
		bool compute_symmetric_gen_partial(void);

		//data
		bool m_full;
		bool m_symmetric;
		uint32_t m_order;
		uint32_t m_index_min;
		uint32_t m_index_max;
		double *m_A, *m_B, *m_sr, *m_si, *m_U, *m_V;
	};
}