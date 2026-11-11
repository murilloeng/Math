#pragma once

//std
#include <cstdint>

//Math
#include "Math/inc/Solvers/Implicit.hpp"
#include "Math/inc/Solvers/Incremental.hpp"

namespace math
{
	namespace solvers
	{
		class Newmark : virtual public Implicit, virtual public Incremental
		{
		public:
			//constructors
			Newmark(void);

			//destructor
			~Newmark(void) = default;

			//solve
			void compute_step(void) override;

			//data
			uint32_t states_set(void) const override;
			uint32_t energies_set(void) const override;
			uint32_t tangents_set(void) const override;

			//types
			typedef std::function<void(double&, double&, double&, const double*, const double*, double)> EnergyCallback;
			typedef std::function<void(double*, double*, double*, double*, double*, const double*, const double*, const double*, double)> SystemCallback;

			//data
			double beta(double);
			double beta(void) const;

			double gamma(double);
			double gamma(void) const;

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
			double m_g, m_b;
			SystemCallback m_system;
			EnergyCallback m_energy;
		};
	}
}