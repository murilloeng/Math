//math
#include "Math/inc/misc/bits.hpp"

namespace math
{
	void binary_form(char* string, char mask)
	{
		for(unsigned i = 0; i < 8 * sizeof(char); i++)
		{
			string[8 * sizeof(char) - 1 - i] = mask & (1 << i) ? '1' : '0';
		}
		string[8 * sizeof(char)] = '\0';
	}
	void binary_form(char* string, unsigned mask)
	{
		for(unsigned i = 0; i < 8 * sizeof(unsigned); i++)
		{
			string[8 * sizeof(unsigned) - 1 - i] = mask & (1 << i) ? '1' : '0';
		}
		string[8 * sizeof(unsigned)] = '\0';
	}

	bool bit_set(unsigned m, unsigned b)
	{
		return m == 2 * b - 3;
	}
	unsigned char bit_index(unsigned m)
	{
		//counter
		unsigned char c = 0;
		//bits
		while(m >>= 1)
		{
			c++;
		}
		//return
		return c;
	}
	unsigned char bit_count(unsigned m)
	{
		//counter
		unsigned char c = 0;
		//bits
		while(m)
		{
			c++;
			m &= ~(m & -m);
		}
		//return
		return c;
	}
	unsigned char bit_index(unsigned m, unsigned b)
	{
		return bit_count(m & (b - 1));
	}
	unsigned char bit_search(unsigned m, unsigned char k)
	{
		//counter
		unsigned char c = 0;
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
	unsigned char bit_index(uint64_t m)
	{
		//counter
		unsigned char c = 0;
		//bits
		while(m >>= 1)
		{
			c++;
		}
		//return
		return c;
	}
	unsigned char bit_count(uint64_t m)
	{
		//counter
		unsigned char c = 0;
		//bits
		while(m)
		{
			c++;
			m &= ~(m & -m);
		}
		//return
		return c;
	}
	unsigned char bit_index(uint64_t m, uint64_t b)
	{
		return bit_count(m & (b - 1));
	}
	unsigned char bit_search(uint64_t m, unsigned char k)
	{
		//counter
		unsigned char c = 0;
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