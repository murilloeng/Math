//std
#include <cmath>

//Math
#include "Math/Test/inc/validation.hpp"
#include "Math/inc/Miscellaneous/util.hpp"
#include "Math/inc/Validation/Validator.hpp"

void tests::validation::validation_from_function(void)
{
	//data
	srand(time(nullptr));
	const uint32_t nn = 100;
	math::validation::Validator validator;
	//Setup
	double xn[nn], yn[nn];
	for(uint32_t i = 0; i < nn; i++)
	{
		const double s = math::randu();
		xn[i] = 2 * M_PI * (i + s) / nn;
		yn[i] = sin(xn[i]);
	}
	//items
	validator.create_item();
	validator.item(0)->function(sin);
	validator.item(0)->load_numeric(xn, yn, nn);
	//validation
	validator.validate();
}