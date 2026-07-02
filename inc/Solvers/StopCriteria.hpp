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
				StepMaximum			= 1 << 0,
				TimeMaximum			= 1 << 1,
				LoadLimitMinimum	= 1 << 2,
				LoadLimitMaximum	= 1 << 3,
				LoadLocalMinimum	= 1 << 4,
				LoadLocalMaximum	= 1 << 5,
				LoadValueNegative	= 1 << 6,
				LoadValuePositive	= 1 << 7,
				StateLimitMinimum	= 1 << 8,
				StateLimitMaximum	= 1 << 9,
				StateLocalMinimum	= 1 << 10,
				StateLocalMaximum	= 1 << 11,
				StateValueNegative	= 1 << 12,
				StateValuePositive	= 1 << 13,
				Last
			};

			//constructor
			StopCriteria(Incremental*);

			//destructor
			~StopCriteria(void);

			//data
			uint32_t types(Type);
			uint32_t types(uint32_t);
			uint32_t types(void) const;

			double load_min(double);
			double load_min(void) const;

			double load_max(double);
			double load_max(void) const;

			double state_min(double);
			double state_min(void) const;

			double state_max(double);
			double state_max(void) const;

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

			//data
			Type m_stop;
			uint32_t m_types;
			Incremental* m_solver;
			double m_p_min, m_p_max;
			double m_x_min, m_x_max;
		};
	}
}