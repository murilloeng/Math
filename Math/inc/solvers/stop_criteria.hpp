#pragma once

//std
#include <cstdint>

namespace math
{
	namespace solvers
	{
		class solver;
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
				step_maximum			= 1 << 0,
				time_maximum			= 1 << 1,
				load_limit_minimum		= 1 << 2,
				load_limit_maximum		= 1 << 3,
				load_local_minimum		= 1 << 4,
				load_local_maximum		= 1 << 5,
				load_value_negative		= 1 << 6,
				load_value_positive		= 1 << 7,
				state_limit_minimum		= 1 << 8,
				state_limit_maximum		= 1 << 9,
				state_local_minimum		= 1 << 10,
				state_local_maximum		= 1 << 11,
				state_value_negative	= 1 << 12,
				state_value_positive	= 1 << 13,
				last
			};

			//constructor
			stop_criteria(void);

			//destructor
			~stop_criteria(void);

			//stop
			bool stop(void);

		private:
			//stop
			bool stop_step_maximum(void) const;
			bool stop_time_maximum(void) const;
			bool stop_load_limit_minimum(void) const;
			bool stop_load_limit_maximum(void) const;
			bool stop_load_local_minimum(void) const;
			bool stop_load_local_maximum(void) const;
			bool stop_load_value_negative(void) const;
			bool stop_load_value_positive(void) const;
			bool stop_state_limit_minimum(void) const;
			bool stop_state_limit_maximum(void) const;
			bool stop_state_local_minimum(void) const;
			bool stop_state_local_maximum(void) const;
			bool stop_state_value_negative(void) const;
			bool stop_state_value_positive(void) const;

		public:
			//data
			type m_stop;
			uint32_t m_types;
			solver* m_solver;
			double m_p_min, m_p_max;
			double m_x_min, m_x_max;
		};
	}
}