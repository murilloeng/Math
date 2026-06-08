#pragma once

//std
#include <cstdint>

namespace math
{
	class Eigen
	{
	public:
		//constructor
		Eigen(void);

		//destructor
		~Eigen(void);

		//enums
		enum class Type : uint32_t
		{
			Full, Index, Value
		};

		//data
		Type type(Type);
		Type type(void) const;

		bool symmetry(bool);
		bool symmetry(void) const;

		uint32_t order(uint32_t);
		uint32_t order(void) const;

		double value_min(double);
		double value_min(void) const;

		double value_max(double);
		double value_max(void) const;

		uint32_t index_min(uint32_t);
		uint32_t index_min(void) const;

		uint32_t index_max(uint32_t);
		uint32_t index_max(void) const;

		const double* data(uint32_t) const;
		const double* data(uint32_t, const double*);

		//compute
		bool compute(bool);

		//eigenvalues
		const double* eigenvalues(uint32_t) const;
		double eigenvalue(uint32_t, uint32_t) const;

		//eigenvectors
		const double* eigenvectors(uint32_t) const;
		const double* eigenvector(uint32_t, uint32_t) const;

	private:
		//setup
		void cleanup(void);
		void allocate(void);

		//compute
		bool compute_symmetric_std_full(bool);
		bool compute_symmetric_std_partial(bool);

		//data
		Type m_type;
		bool m_symmetry;
		uint32_t m_order;
		uint32_t m_modes;
		double m_value_min;
		double m_value_max;
		uint32_t m_index_min;
		uint32_t m_index_max;
		const double* m_data[2];
		double* m_eigenvalues[2];
		double* m_eigenvectors[2];

	};
}