#pragma once

//Math
#include "Math/inc/Groups/LieGroup.hpp"

namespace math
{
	class ASE3;
}

namespace math
{
	class GSE3 : public LieGroup<GSE3, ASE3>
	{
	public:
		//constructor
		explicit GSE3(double, double, double, double, double, double, double);
		explicit GSE3(double = 0, double = 0, double = 0, double = 0, double = 0, double = 0);

		//destructor
		~GSE3(void) = default;

		//random
		GSE3& randu(double, double) override;

		//actions
		Matrix matrix(void) const override;
		Vector action(const Vector&) const override;

		//operations
		ASE3 log(void) const override;
		GSE3 inverse(void) const override;
		GSE3 operator*(const GSE3&) const override;
		GSE3 interpolate(const GSE3&, double) const override;

		//tests
		static bool test_matrix(uint32_t, double);
		static bool test_inverse(uint32_t, double);
		static bool test_log_exp(uint32_t, double);

	private:
		//data
		double m_data[7];
	};
}