//std
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <cstdint>
#include <cstring>

//Math
#include "Math/inc/Linear/Vec3.hpp"
#include "Math/inc/Miscellaneous/BoundingBox.hpp"

namespace math
{
	//constructor
	BoundingBox::BoundingBox(void) : m_min{-DBL_MAX, -DBL_MAX, -DBL_MAX}, m_max{+DBL_MAX, +DBL_MAX, +DBL_MAX}
	{
		return;
	}

	//print
	void BoundingBox::print(const char* header) const
	{
		if(strlen(header)) printf("%s\n", header);
		printf("Min: %+.2e %+.2e %+.2e\n", m_min[0], m_min[1], m_min[2]);
		printf("Max: %+.2e %+.2e %+.2e\n", m_max[0], m_max[1], m_max[2]);
	}

	//data
	double BoundingBox::radius(void) const
	{
		//data
		return sqrt(
			(m_max[0] - m_min[0]) * (m_max[0] - m_min[0]) +
			(m_max[1] - m_min[1]) * (m_max[1] - m_min[1]) +
			(m_max[2] - m_min[2]) * (m_max[2] - m_min[2]));
	}
	const double* BoundingBox::min(void) const
	{
		return m_min;
	}
	const double* BoundingBox::max(void) const
	{
		return m_max;
	}

	//status
	bool BoundingBox::status(void) const
	{
		bool test = true;
		for(const double& x : m_min) test = test && x != -DBL_MAX;
		for(const double& x : m_max) test = test && x != +DBL_MAX;
		return test;
	}

	//compute
	void BoundingBox::reset(void)
	{
		for(double& x : m_min) x = -DBL_MAX;
		for(double& x : m_max) x = +FLT_MAX;
	}
	void BoundingBox::repair(void)
	{
		if(!status())
		{
			for(double& x : m_min) x = -1;
			for(double& x : m_max) x = +1;
		}
		if((Vec3(m_max) - Vec3(m_min)).norm() == 0)
		{
			for(double& x : m_min) x -= 1;
			for(double& x : m_max) x += 1;
		}
	}
	void BoundingBox::insert_vertex(const double* vertex)
	{
		for(uint32_t i = 0; i < 3; i++)
		{
			if(m_min[i] == -DBL_MAX || m_min[i] > vertex[i]) m_min[i] = vertex[i];
			if(m_max[i] == +DBL_MAX || m_max[i] < vertex[i]) m_max[i] = vertex[i];
		}
	}

	//combine
	void BoundingBox::combine(const BoundingBox& box)
	{
		for(uint32_t i = 0; i < 3; i++)
		{
			m_min[i] = fminf(m_min[i], box.m_min[i]);
			m_max[i] = fmaxf(m_max[i], box.m_max[i]);
		}
	}
	BoundingBox BoundingBox::combine(const BoundingBox& box_1, const BoundingBox& box_2)
	{
		//combine
		BoundingBox box;
		for(uint32_t i = 0; i < 3; i++)
		{
			box.m_min[i] = fminf(box_1.m_min[i], box_2.m_min[i]);
			box.m_max[i] = fmaxf(box_1.m_max[i], box_2.m_max[i]);
		}
		//return
		return box;
	}
}