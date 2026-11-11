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
			uint32_t states_set(void) const override;
			uint32_t forces_set(void) const override;
			uint32_t energies_set(void) const override;
			uint32_t tangents_set(void) const override;

			//types
			typedef std::function<void(double&, double&, const double*, double)> EnergyCallback;
			typedef std::function<void(double*, double*, double*, const double*, double)> SystemCallback;

			//data
			EnergyCallback energy(void) const;
			EnergyCallback energy(EnergyCallback);

			SystemCallback system(void) const;
			SystemCallback system(SystemCallback);

		protected:
			//solve
			void check(void) override;
			void print(void) override;
			void setup(void) override;
			void compute(void) override;
			bool predictor(void) override;
			bool corrector(void) override;

			//data
			EnergyCallback m_energy;
			SystemCallback m_system;
		};
	}
}