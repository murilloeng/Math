#pragma once

//std
#include <cstdint>
#include <functional>

//Math
#include "Math/inc/Solvers/Implicit.hpp"
#include "Math/inc/Solvers/Incremental.hpp"

namespace math
{
	namespace solvers
	{
		class NewtonRaphson : virtual public Implicit, virtual public Incremental
		{
		public:
			//constructors
			NewtonRaphson(void);

			//destructor
			~NewtonRaphson(void);

			//solve
			void step(void) override;

			//data
			uint32_t state_set(void) const override;
			uint32_t force_set(void) const override;
			uint32_t tangent_set(void) const override;

			//types
			typedef std::function<void(double*, double, const double*)> Residue;
			typedef std::function<void(double*, double, const double*)> Tangent_p;
			typedef std::function<void(double*, double, const double*)> Tangent_x;
			typedef std::function<void(double*, double*, const double*)> System_1;
			typedef std::function<void(double*, double*, double*, double, const double*)> System_2;

			//data
			Residue residue(Residue);
			Residue residue(void) const;

			System_1 system_1(System_1);
			System_1 system_1(void) const;

			System_2 system_2(System_2);
			System_2 system_2(void) const;

			Tangent_p tangent_p(Tangent_p);
			Tangent_p tangent_p(void) const;

			Tangent_x tangent_x(Tangent_x);
			Tangent_x tangent_x(void) const;

		protected:
			//solve
			void check(void) override;
			void print(void) override;
			void setup(void) override;
			void compute(void) override;
			void predictor(void) override;
			void corrector(void) override;

			//data
			Residue m_residue;
			System_1 m_system_1;
			System_2 m_system_2;
			Tangent_p m_tangent_p;
			Tangent_x m_tangent_x;
		};
	}
}