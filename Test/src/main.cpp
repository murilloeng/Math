//std
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdint>
#include <cstdlib>

//math
#include "Math/Math/inc/linear/vector.hpp"
#include "Math/Math/inc/geometry/Circle.hpp"

//test
#include "Math/Test/inc/tests.hpp"

int main(void)
{
	//test
	// tests::solvers::harmonic::pyramid();

	double r, q1, q2, q3;
	math::vec3 x1, x2, x3, xc, n, t1, t2;
	
	srand(uint32_t(time(nullptr)));
	for (uint32_t i = 0; i < 1000; i++)
	{
		n.randu();
		xc.randu();
		n.normalize();
		math::vector(&r, 1).randu(0, 1);
		math::vector(&q1, 1).randu(0, 1);
		math::vector(&q2, 1).randu(0, 1);
		math::vector(&q3, 1).randu(0, 1);
		n.triad(t1, t2);

		x1 = xc + r * (cos(q1) * t1 + sin(q1) * t2);
		x2 = xc + r * (cos(q2) * t1 + sin(q2) * t2);
		x3 = xc + r * (cos(q3) * t1 + sin(q3) * t2);
		math::geometry::Circle circle(x1, x2, x3);

		if(fabs(r - circle.m_radius) > 1e-5)
		{
			printf("Error: radius!\n");
			break;
		}
		if((n - circle.m_normal).norm() > 1e-5 && (n + circle.m_normal).norm() > 1e-5)
		{
			printf("Error: normal!\n");
			break;
		}
		if((xc - circle.m_center).norm() > 1e-5)
		{
			printf("Error: center!\n");
			break;
		}
		printf("%d: ok\n", i);
	}

	//return
	return EXIT_SUCCESS;
}