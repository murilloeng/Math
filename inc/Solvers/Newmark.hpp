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
			~Newmark(void);

			//solve
			void compute_step(void) override;

			//data
			uint32_t state_set(void) const override;
			uint32_t force_set(void) const override;
			uint32_t tangent_set(void) const override;

			//types
			typedef std::function<void(double*, const double*, const double*)> InternalForce;
			typedef std::function<void(double*, const double*, const double*, double)> ExternalForce;

			typedef std::function<void(double*, const double*)> Inertia;
			typedef std::function<void(double*, const double*, const double*, double)> Damping;
			typedef std::function<void(double*, const double*, const double*, const double*, double)> Stiffness;

			//data
			double beta(double);
			double beta(void) const;

			double gamma(double);
			double gamma(void) const;

			Inertia inertia(Inertia);
			Inertia inertia(void) const;

			Damping damping(Damping);
			Damping damping(void) const;

			Stiffness stiffness(Stiffness);
			Stiffness stiffness(void) const;

			InternalForce internal_force(void) const;
			InternalForce internal_force(InternalForce);

			ExternalForce external_force(void) const;
			ExternalForce external_force(ExternalForce);

		private:
			//solve
			void check(void) override;
			void print(void) override;
			void setup(void) override;
			void compute(void) override;
			void predictor(void) override;
			void corrector(void) override;

			//data
			double m_g, m_b;

			Inertia m_inertia;
			Damping m_damping;
			Stiffness m_stiffness;

			InternalForce m_internal_force;
			ExternalForce m_external_force;
		};
	}
}