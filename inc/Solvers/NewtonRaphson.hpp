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
			NewtonRaphson(void) = default;

			//destructor
			~NewtonRaphson(void) = default;

			//solve
			void compute_step(void) override;

			//data
			uint32_t state_set(void) const override;
			uint32_t force_set(void) const override;
			uint32_t energy_set(void) const override;
			uint32_t tangent_set(void) const override;

			//types
			typedef std::function<void(double&, double&, double, const double*)> Energies;
			typedef std::function<void(double*, double*, double*, double, const double*)> System;

			//data
			System system(System);
			System system(void) const;

			Energies energies(Energies);
			Energies energies(void) const;

		protected:
			//solve
			void check(void) override;
			void print(void) override;
			void setup(void) override;
			void compute(void) override;
			bool predictor(void) override;
			bool corrector(void) override;

			//data
			System m_system;
			Energies m_energies;
		};
	}
}