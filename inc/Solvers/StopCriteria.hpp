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
				StateLimitMinimum,
				StateLimitMaximum,
				StateLocalMinimum,
				StateLocalMaximum,
				StateValueNegative,
				StateValuePositive,
				ParameterLimitMinimum,
				ParameterLimitMaximum,
				ParameterLocalMinimum,
				ParameterLocalMaximum,
				ParameterValueNegative,
				ParameterValuePositive,
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

			double state_min(double);
			double state_min(void) const;

			double state_max(double);
			double state_max(void) const;

			double parameter_min(double);
			double parameter_min(void) const;

			double parameter_max(double);
			double parameter_max(void) const;

		private:
			//stop
			bool stop_step_maximum(void) const;
			bool stop_time_maximum(void) const;
			bool stop_state_limit_minimum(void) const;
			bool stop_state_limit_maximum(void) const;
			bool stop_state_local_minimum(void) const;
			bool stop_state_local_maximum(void) const;
			bool stop_state_value_negative(void) const;
			bool stop_state_value_positive(void) const;
			bool stop_parameter_limit_minimum(void) const;
			bool stop_parameter_limit_maximum(void) const;
			bool stop_parameter_local_minimum(void) const;
			bool stop_parameter_local_maximum(void) const;
			bool stop_parameter_value_negative(void) const;
			bool stop_parameter_value_positive(void) const;

			//data
			Type m_stop;
			uint32_t m_types;
			Incremental* m_solver;
			double m_state_min, m_parameter_min;
			double m_state_max, m_parameter_max;
		};
	}
}