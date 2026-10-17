#pragma once

//Math
#include "Math/inc/Groups/LieGroup.hpp"

namespace math
{
	class ASO2;
}

namespace math
{
	class GSO2 : public LieGroup<GSO2, ASO2>
	{
	public:
		//constructor
		explicit GSO2(double = 0);

		//destructor
		~GSO2(void) = default;

		//actions
		Matrix matrix(void) const override;
		Vector action(const Vector&) const override;

		//operations
		ASO2 log(void) const override;
		GSO2 inverse(void) const override;
		GSO2 operator*(const GSO2&) const override;
		GSO2 interpolate(const GSO2&, double) const override;

	private:
		//data
		double m_data;
	};
}