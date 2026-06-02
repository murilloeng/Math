//std
#include <cstring>

//Math
#include "Math/inc/Validation/Validator.hpp"

namespace math
{
	namespace validation
	{
		//constructor
		Validator::Validator(void) : m_silent{false}
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

		//data
		bool Validator::silent(void) const
		{
			return m_silent;
		}
		bool Validator::silent(bool silent)
		{
			for(Item* item : m_items)
			{
				item->m_silent = silent;
			}
			return m_silent = silent;
		}

		void Validator::create_item(void)
		{
			m_items.push_back(new Item);
		}
		Item* Validator::item(uint32_t index) const
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
			if(!m_silent) printf("Status: %s\n", test ? "✅" : "❌");
			//return
			return test;
		}
	}
}