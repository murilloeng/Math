#pragma once

//Math
#include "Math/inc/Solvers/Incremental.hpp"

namespace math
{
	namespace solvers
	{
		class RungeKutta : virtual public Incremental
		{
		public:
			//constructors
			RungeKutta(void) = default;

			//destructor
			~RungeKutta(void) = default;

			//solve
			void compute_step(void) override;

			//data
			uint32_t states_set(void) const override;
			uint32_t energies_set(void) const override;
			uint32_t tangents_set(void) const override;

			//types
			typedef std::function<void(double&, double&, double&, const double*, const double*, double)> EnergyCallback;
			typedef std::function<void(double*, double*, double*, const double*, const double*, double)> SystemCallback;

			//data
			EnergyCallback energy(void) const;
			EnergyCallback energy(EnergyCallback);

			SystemCallback system(void) const;
			SystemCallback system(SystemCallback);

		private:
			//solve
			void check(void) override;
			void print(void) override;
			void setup(void) override;
			void compute(void) override;
			bool predictor(void) override;
			bool corrector(void) override;

			//compute
			void compute_tangent_1(void);
			void compute_tangent_2(void);
			void compute_tangent_3(void);
			void compute_tangent_4(void);

			//data
			EnergyCallback m_energy;
			SystemCallback m_system;
		};
	}
}