#pragma once

//std
#include <cstdint>

namespace math
{
	namespace solvers
	{
		class newton_raphson;
	}
}

namespace math
{
	namespace solvers
	{
		class stop_criteria
		{
		public:
			//types
			enum class type : uint32_t
			{
				load_maximum	= 1 << 0,
				load_minimum	= 1 << 1,
				load_positive	= 1 << 2,
				load_negative	= 1 << 3,
				state_maximum	= 1 << 4,
				state_minimum	= 1 << 5,
				state_positive	= 1 << 6,
				state_negative	= 1 << 7,
				last
			};

			//constructor
			stop_criteria(void);

			//destructor
			~stop_criteria(void);

			//stop
			bool stop(void) const;

		private:
			//stop
			bool stop_load_maximum(void) const;
			bool stop_load_minimum(void) const;
			bool stop_load_positive(void) const;
			bool stop_load_negative(void) const;
			bool stop_state_maximum(void) const;
			bool stop_state_minimum(void) const;
			bool stop_state_positive(void) const;
			bool stop_state_negative(void) const;

		public:
			//data
			uint32_t m_types;
			double m_p_min, m_p_max;
			double m_x_min, m_x_max;
			newton_raphson* m_solver;
		};
	}
}