//std
#include <cstring>

//Math
#include "Math/inc/Validation/Validator.hpp"

namespace math
{
	namespace validation
	{
		//constructor
		Validator::Validator(void)
		{
			return;
		}

		//destructor
		Validator::~Validator(void)
		{
			for(const Item* item : m_items)
			{
				delete item;
			}
		}

		//items
		const Item* Validator::item(uint32_t index) const
		{
			return m_items[index];
		}

		//validation
		bool Validator::validate(void)
		{
			//validation
			bool test = true;
			for(Item* item : m_items)
			{
				test = test && item->validate();
			}
			printf("Status: %s\n", test ? "✅" : "❌");
			//return
			return test;
		}
	}
}