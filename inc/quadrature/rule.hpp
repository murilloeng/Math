#pragma once

//std
#include <cstdint>

namespace math
{
	namespace quadrature
	{
		enum class rule : uint32_t
		{
			Lobatto,
			Legendre,
			Last
		};
	}
}