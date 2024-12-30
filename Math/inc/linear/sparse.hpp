#pragma once

namespace math
{
	class sparse
	{
	public:
		//constructors
		sparse(uint32_t);

		//destructor
		~sparse(void);

	private:
		//data
		double* m_data;
		uint32_t m_size;
		uint32_t* m_col_map;
		uint32_t* m_row_map;
	};
}