//std
#include <ctime>
#include <cmath>
#include <cstdio>

//Math
#include "Math/inc/Linear/Vec3.hpp"
#include "Math/inc/Linear/Mat3.hpp"

//Test
#include "Math/Test/inc/rotations.hpp"

void tests::rotations::Vec3::rotation_tensor(void)
{
	//data
	const uint32_t nt = 10000;
	math::Vec3 v, s[3], r1[3], r2[3];
	//test
	srand((uint32_t) time(nullptr));
	for(uint32_t i = 0; i < nt; i++)
	{
		v.randu();
		s[0].randu();
		s[0].normalize();
		s[0].triad(s[1], s[2]);
		for(uint32_t j = 0; j < 3; j++)
		{
			r1[j] = v.rotate(s[j]);
			r2[j] = v.rotation_tensor() * s[j];
		}
		bool test = true;
		for(uint32_t j = 0; j < 3; j++)
		{
			const uint32_t k = (j + 1) % 3;
			test = test && fabs(r1[j].inner(r1[k])) < 1e-5;
			test = test && fabs(r2[j].inner(r2[k])) < 1e-5;
			test = test && fabs(r1[j].inner(r1[j]) - 1) < 1e-5;
			test = test && fabs(r2[j].inner(r2[j]) - 1) < 1e-5;
		}
		printf("Test %04d: %s\n", i, test ? "ok" : "not ok");
		if(!test) break;
	}
}