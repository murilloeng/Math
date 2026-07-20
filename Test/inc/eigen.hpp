#pragma once

namespace tests
{
	namespace eigen
	{
		void dense_svd(void);
		void dense_sym_std_full(void);
		void dense_sym_gen_full(void);
		void dense_sym_std_part(void);
		void dense_sym_gen_part(void);
		void dense_non_std_full(void);
		void dense_non_gen_full(void);

		void sparse_svd(void);
		void sparse_sym_std_part(void);
		void sparse_sym_gen_part(void);
		void sparse_non_std_part(void);
		void sparse_non_gen_part(void);
	}
}