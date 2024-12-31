//math
#include "Galileo/mat/inc/misc/bits.hpp"

namespace math
{
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
	uint32_t bit_search(uint32_t m, uint8_t k)
	{
		//bits
		for(uint8_t c = 0; c < k; c++)
		{
			m &= ~(m & -m);
		}
		//return
		return m & -m;
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
	uint64_t bit_search(uint64_t m, uint8_t k)
	{
		//bits
		for(uint8_t c = 0; c < k; c++)
		{
			m &= ~(m & -m);
		}
		//return
		return m & -m;
	}
}