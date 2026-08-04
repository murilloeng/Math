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

			//data
			uint32_t add_type(Type);
			uint32_t types(void) const;

			Type stop_type(void) const;
			
			double load_min(double);
			double load_min(void) const;
			
			double load_max(double);
			double load_max(void) const;

			double state_min(double);
			double state_min(void) const;

			double state_max(double);
			double state_max(void) const;

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
			Type m_stop_type;
			uint32_t m_types;
			Incremental* m_solver;
			double m_state_min, m_load_min;
			double m_state_max, m_load_max;
		};
	}
}