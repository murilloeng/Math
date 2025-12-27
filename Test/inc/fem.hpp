#pragma once

namespace tests
{
	namespace fem
	{
		namespace beam
		{
			namespace dynamics
			{
				void strains(void);
				void rotation(void);
			}
		}
		void beam3DCR(void);
		void beam3DTL(void);
		void revolute_fixed(void);
		void revolute_flexible(void);
	}
}