#pragma once

//std
#include <cstdint>

namespace math
{
	namespace solvers
	{
		class Incremental;
	}
}

namespace math
{
	namespace solvers
	{
		class StopCriteria
		{
		public:
			//types
			enum class Type : uint32_t
			{
				StepMaximum,
				TimeMaximum,
				LoadLimitMinimum,
				LoadLimitMaximum,
				LoadLocalMinimum,
				LoadLocalMaximum,
				LoadValueNegative,
				LoadValuePositive,
				StateLimitMinimum,
				StateLimitMaximum,
				StateLocalMinimum,
				StateLocalMaximum,
				StateValueNegative,
				StateValuePositive,
				Last
			};

			//constructor
			StopCriteria(Incremental*);

			//destructor
			~StopCriteria(void);

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
			Type m_stop;
			uint32_t m_types;
			Incremental* m_solver;
			double m_p_min, m_p_max;
			double m_x_min, m_x_max;
		};
	}
}