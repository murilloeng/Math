#pragma once

namespace tests
{
	namespace eigen
	{
		void dense_symmetric_std_full(void);
		void dense_symmetric_gen_full(void);
		void dense_symmetric_std_partial(void);
		void dense_symmetric_gen_partial(void);
		void dense_non_symmetric_std_full(void);
		void dense_non_symmetric_gen_full(void);
		void dense_singular_value_decomposition(void);
		
		void sparse_symmetric_std_partial(void);
		void sparse_symmetric_gen_partial(void);
		void sparse_non_symmetric_std_partial(void);
		void sparse_non_symmetric_gen_partial(void);
		void sparse_singular_value_decomposition(void);
	}
}