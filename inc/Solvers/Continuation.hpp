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
				MinimalNorm,
				LoadControl,
				StateControl,
				ArcLengthSpherical,
				ArcLengthCylindrical,
				Last
			};

			//constructor
			Continuation(Solver*, Type);

			//destructor
			~Continuation(void);

			//continuation
			double predictor(void) const;
			double corrector(void) const;

		private:
			//continuation
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

		public:
			//data
			Type m_type;
			Solver* m_solver;
		};
	}
}