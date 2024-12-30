#pragma once

//std
#include <cstdint>

namespace math
{
	enum class rule : uint32_t
	{
		lobatto,
		legendre,
		last
	};
}