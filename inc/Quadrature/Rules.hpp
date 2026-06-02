#pragma once

//std
#include <cstdint>

namespace math
{
	namespace quadrature
	{
		enum class Rule : uint32_t
		{
			Lobatto,
			Legendre,
			Last
		};
	}
}