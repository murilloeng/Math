#pragma once

namespace math
{
	enum class strategy : uint32_t
	{
		arc_length		= 1 << 0,
		minimal_norm	= 1 << 1,
		control_load	= 1 << 2,
		control_state	= 1 << 3,
		last
	};
}