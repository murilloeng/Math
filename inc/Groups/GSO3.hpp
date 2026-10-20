#pragma once

//Math
#include "Math/inc/Groups/LieGroup.hpp"

namespace math
{
	class ASO3;
}

namespace math
{
	class GSO3 : public LieGroup<GSO3, ASO3>
	{
	public:
		//constructor
		explicit GSO3(const double*);
		explicit GSO3(double = 1, double = 0, double = 0, double = 0);

		//destructor
		~GSO3(void) = default;

		//random
		GSO3& randu(double, double) override;

		//actions
		Matrix matrix(void) const override;
		Vector action(const Vector&) const override;

		//operations
		ASO3 log(void) const override;
		GSO3 inverse(void) const override;
		GSO3 operator*(const GSO3&) const override;
		GSO3 interpolate(const GSO3&, double) const override;

		double& operator[](uint32_t);
		const double& operator[](uint32_t) const;

		//tests
		static bool test_matrix(uint32_t, double);
		static bool test_inverse(uint32_t, double);
		static bool test_log_exp(uint32_t, double);

	private:
		//data
		double m_data[4];
	};
}