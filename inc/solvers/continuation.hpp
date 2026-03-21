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
		class continuation
		{
		public:
			//types
			enum class type : uint32_t
			{
				minimal_norm			= 1 << 0,
				control_load			= 1 << 1,
				control_state			= 1 << 2,
				arc_length_spherical	= 1 << 3,
				arc_length_cylindrical	= 1 << 4,
				last
			};
	
			//constructor
			continuation(void);
			continuation(type);
	
			//destructor
			~continuation(void);
	
			//continuation
			double predictor(void) const;
			double corrector(void) const;
	
			//types
			double predictor_minimal_norm(void) const;
			double corrector_minimal_norm(void) const;
			double predictor_control_load(void) const;
			double corrector_control_load(void) const;
			double predictor_control_state(void) const;
			double corrector_control_state(void) const;
			double predictor_arc_length_spherical(void) const;
			double corrector_arc_length_spherical(void) const;
			double predictor_arc_length_cylindrical(void) const;
			double corrector_arc_length_cylindrical(void) const;
	
			//data
			type m_type;
			solver* m_solver;
		};
	}
}