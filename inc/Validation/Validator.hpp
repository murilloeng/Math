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
	
			//types
			typedef std::vector<uint32_t> ulist;
	
			//items
			const Item* item(uint32_t) const;
	
			//validation
			bool validate(void);
	
		private:
			//data
			std::vector<Item*> m_items;
		};
	}
}