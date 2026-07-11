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
			RungeKutta(void);

			//destructor
			~RungeKutta(void);

			//solve
			void compute_step(void) override;

			//data
			uint32_t state_set(void) const override;
			uint32_t force_set(void) const override;
			uint32_t tangent_set(void) const override;

			//types
			typedef std::function<void(double*, const double*)> Inertia;
			typedef std::function<void(double*, const double*, const double*)> InternalForce;
			typedef std::function<void(double*, const double*, const double*, double)> ExternalForce;

			//data
			Inertia inertia(Inertia);
			Inertia inertia(void) const;

			InternalForce internal_force(void) const;
			InternalForce internal_force(InternalForce);

			ExternalForce external_force(void) const;
			ExternalForce external_force(ExternalForce);

		private:
			//solve
			void check(void) override;
			void print(void) override;
			void compute(void) override;
			void predictor(void) override;
			void corrector(void) override;

			//compute
			void compute_tangent_1(void);
			void compute_tangent_2(void);
			void compute_tangent_3(void);
			void compute_tangent_4(void);

			//data
			Inertia m_inertia;
			InternalForce m_internal_force;
			ExternalForce m_external_force;

		};
	}
}