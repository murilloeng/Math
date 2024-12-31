#pragma once

//std
#include <cstdint>

namespace math
{
	namespace solvers
	{
		class newmark
		{
		public:
			//constructors
			newmark(uint32_t, bool);

			//destructor
			~newmark(void);

		private:
			//data
			void update(void);
			void residue(void);
			void inertia(void);
			void damping(void);
			void stifness(void);
			void internal(void);
			void external(void);

			//solve
			void setup(void);
			void predictor(void);
			void corrector(void);
			void serialize(void);

		public:
			//solve
			void step(void);
			void solve(void);

			//data
			bool m_mem;
			uint32_t m_s, m_nd, m_ns;
			double m_t, m_T, m_g, m_b, *m_x, *m_v, *m_a, *m_r;
			double m_dt, *m_dx, *m_dv, *m_fe, *m_fi, *m_K, *m_C, *m_M;

			void (*m_external)(double*, double);
			void (*m_internal)(double*, const double*, const double*);

			void (*m_inertia)(double*, const double*);
			void (*m_damping)(double*, const double*, const double*);
			void (*m_stiffness)(double*, const double*, const double*);
		};
	}
}