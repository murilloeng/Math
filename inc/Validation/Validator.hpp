#pragma once

//std
#include <vector>
#include <string>

//Math
#include "Math/inc/Validation/Item.hpp"

namespace math
{
	namespace validation
	{
		class Validator
		{
		public:
			//constructor
			Validator(void);
	
			//destructor
			~Validator(void);
	
			//items
			void create_item(void);
			Item* item(uint32_t) const;
	
			//validation
			bool validate(void);
	
		private:
			//data
			std::vector<Item*> m_items;
		};
	}
}