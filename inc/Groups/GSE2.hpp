#pragma once

//Math
#include "Math/inc/Groups/LieGroup.hpp"

namespace math
{
	class ASE2;
}

namespace math
{
	class GSE2 : public LieGroup<GSE2, ASE2>
	{
	public:
		//constructor
		explicit GSE2(const double*);
		explicit GSE2(double = 0, double = 0, double = 0);

		//destructor
		~GSE2(void) = default;

		//random
		GSE2& randu(double, double) override;

		//actions
		Matrix matrix(void) const override;
		Vector action(const Vector&) const override;

		//operators
		double& operator[](uint32_t) override;
		const double& operator[](uint32_t) const override;

		//group
		ASE2 log(void) const override;
		GSE2 inverse(void) const override;
		GSE2 operator*(const GSE2&) const override;
		GSE2 interpolate(const GSE2&, double) const override;

		//tests
		static bool test_matrix(uint32_t, double);
		static bool test_inverse(uint32_t, double);
		static bool test_log_exp(uint32_t, double);

	private:
		//data
		double m_data[3];
	};
}