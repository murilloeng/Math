#pragma once

//std
#include <cstdint>

namespace math
{
	namespace solvers
	{
		class Solver;
	}
}

namespace math
{
	namespace solvers
	{
		class Continuation
		{
		public:
			//types
			enum class Type : uint32_t
			{
				MinimalNorm				= 1 << 0,
				LoadControl				= 1 << 1,
				StateControl			= 1 << 2,
				ArcLengthSpherical		= 1 << 3,
				ArcLengthCylindrical	= 1 << 4,
				Last
			};
	
			//constructor
			Continuation(void);
			Continuation(Type);
	
			//destructor
			~Continuation(void);
	
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
			Type m_type;
			Solver* m_solver;
		};
	}
}