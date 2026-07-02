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
		class Convergence
		{
		public:
			//types
			enum class Type : uint32_t
			{
				Fixed = 1 << 0, 
				Force = 1 << 1,
				Last
			};

			//constructor
			Convergence(Solver*);

			//destructor
			~Convergence(void);

			//data
			Type type(Type);
			Type type(void) const;

			double tolerance(double);
			double tolerance(void) const;

			//check
			bool check(void) const;

		private:
			//data
			Type m_type;
			Solver* m_solver;
			double m_tolerance;
		};
	}
}