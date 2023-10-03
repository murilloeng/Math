#pragma once

namespace math
{
	class sparse
	{
	public:
		//constructors
		sparse(unsigned);

		//destructor
		~sparse(void);

	private:
		//data
		double* m_data;
		unsigned m_size;
		unsigned* m_col_map;
		unsigned* m_row_map;
	};
}