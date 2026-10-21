//std
#include <cstdio>

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
			Item* item = new Item;
			item->m_silent = m_silent;
			m_items.push_back(item);
		}
		Item* Validator::item(uint32_t index) const
		{
			return m_items[index];
		}

		//validation
		bool Validator::validate(void)
		{
			//validation
			for(uint32_t i = 0; i < m_items.size(); i++)
			{
				if(!m_items[i]->validate())
				{
					if(!m_silent) printf("Status: ❌ - Item %u failed!\n", i);
					return false;
				}
			}
			if(!m_silent) printf("Status: ✅\n");
			//return
			return true;
		}
	}
}