#pragma once

//std
#include <cstdint>

namespace math
{
	namespace solvers
	{
		class Implicit;
	}
}

namespace math
{
	namespace solvers
	{
		class Convergence
		{
		public:
			//types
			enum class Type : uint32_t
			{
				Fixed, Force, State, Last
			};

			//constructor
			Convergence(Implicit*);

			//destructor
			~Convergence(void) = default;

			//data
			uint32_t types(uint32_t);
			uint32_t types(void) const;

			uint32_t type_set(Type);
			uint32_t type_add(Type);
			uint32_t type_remove(Type);

			double tolerance(double);
			double tolerance(void) const;

			//check
			bool check(void) const;

		private:
			//data
			uint32_t m_types;
			Implicit* m_solver;
			double m_tolerance;
		};
	}
}