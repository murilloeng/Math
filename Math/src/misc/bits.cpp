//math
#include "Math/Math/inc/misc/bits.hpp"

namespace math
{
	void binary_form(char* string, char mask)
	{
		for(uint32_t i = 0; i < 8 * sizeof(char); i++)
		{
			string[8 * sizeof(char) - 1 - i] = mask & (1 << i) ? '1' : '0';
		}
		string[8 * sizeof(char)] = '\0';
	}
	void binary_form(char* string, uint32_t mask)
	{
		for(uint32_t i = 0; i < 8 * sizeof(uint32_t); i++)
		{
			string[8 * sizeof(uint32_t) - 1 - i] = mask & (1 << i) ? '1' : '0';
		}
		string[8 * sizeof(uint32_t)] = '\0';
	}

	bool bit_set(uint32_t m, uint32_t b)
	{
		return m == 2 * b - 3;
	}
	uint8_t bit_index(uint32_t m)
	{
		//counter
		uint8_t c = 0;
		//bits
		while(m >>= 1)
		{
			c++;
		}
		//return
		return c;
	}
	uint8_t bit_count(uint32_t m)
	{
		//counter
		uint8_t c = 0;
		//bits
		while(m)
		{
			c++;
			m &= ~(m & -m);
		}
		//return
		return c;
	}
	uint8_t bit_index(uint32_t m, uint32_t b)
	{
		return bit_count(m & (b - 1));
	}
	uint8_t bit_search(uint32_t m, uint8_t k)
	{
		//counter
		uint8_t c = 0;
		//bits
		while(c < k)
		{
			c++;
			m &= ~(m & -m);
		}
		//return
		return bit_index(m & -m);
	}

	bool bit_set(uint64_t m, uint64_t b)
	{
		return m == 2 * b - 3;
	}
	uint8_t bit_index(uint64_t m)
	{
		//counter
		uint8_t c = 0;
		//bits
		while(m >>= 1)
		{
			c++;
		}
		//return
		return c;
	}
	uint8_t bit_count(uint64_t m)
	{
		//counter
		uint8_t c = 0;
		//bits
		while(m)
		{
			c++;
			m &= ~(m & -m);
		}
		//return
		return c;
	}
	uint8_t bit_index(uint64_t m, uint64_t b)
	{
		return bit_count(m & (b - 1));
	}
	uint8_t bit_search(uint64_t m, uint8_t k)
	{
		//counter
		uint8_t c = 0;
		//bits
		while(c < k)
		{
			c++;
			m &= ~(m & -m);
		}
		//return
		return bit_index(m & -m);
	}
}