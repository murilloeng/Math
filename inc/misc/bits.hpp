#pragma once

//std
#include <cstdint>

namespace math
{
	void binary_form(char*, char);
	void binary_form(char*, unsigned);

	bool bit_set(unsigned, unsigned);
	unsigned char bit_index(unsigned);
	unsigned char bit_count(unsigned);
	unsigned char bit_index(unsigned, unsigned);
	unsigned char bit_search(unsigned, unsigned char);

	bool bit_set(uint64_t, uint64_t);
	unsigned char bit_index(uint64_t);
	unsigned char bit_count(uint64_t);
	unsigned char bit_index(uint64_t, uint64_t);
	unsigned char bit_search(uint64_t, unsigned char);

}