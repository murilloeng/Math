#pragma once

//Math
#include "Math/inc/Groups/GSO3.hpp"
#include "Math/inc/Groups/LieGroup.hpp"

namespace math
{
	class ASE3;
}

namespace math
{
	class GSE3 : public LieGroup<GSE3, ASE3, 4, 3>
	{
	public:
		//constructor
		explicit GSE3(const double*);
		explicit GSE3(double = 0, double = 0, double = 0, double = 1, double = 0, double = 0, double = 0);

		//destructor
		~GSE3(void) = default;

		//random
		GSE3& randu(double, double) override;

		//actions
		GSO3 gso3(void) const;
		Mat4 matrix(void) const override;
		Vec3 action(const Vec3&) const override;

		//operators
		double& operator[](uint32_t) override;
		const double& operator[](uint32_t) const override;

		//group
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