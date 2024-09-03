#pragma once

//std
#include <cstdint>

namespace math
{
	void binary_form(char*, char);
	void binary_form(char*, uint32_t);

	bool bit_set(uint32_t, uint32_t);
	uint8_t bit_index(uint32_t);
	uint8_t bit_count(uint32_t);
	uint8_t bit_index(uint32_t, uint32_t);
	uint8_t bit_search(uint32_t, uint8_t);

	bool bit_set(uint64_t, uint64_t);
	uint8_t bit_index(uint64_t);
	uint8_t bit_count(uint64_t);
	uint8_t bit_index(uint64_t, uint64_t);
	uint8_t bit_search(uint64_t, uint8_t);

}