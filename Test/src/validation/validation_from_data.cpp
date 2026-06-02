//std
#include <cmath>
#include <ctime>

//Math
#include "Math/Test/inc/validation.hpp"
#include "Math/inc/Miscellaneous/util.hpp"
#include "Math/inc/Validation/Validator.hpp"

void tests::validation::validation_from_data(void)
{
	//data
	srand(time(nullptr));
	const uint32_t nr = 100;
	const uint32_t nn = 100;
	double xn[5 * nn], xr[6 * nr];
	for(uint32_t i = 0; i < nn; i++)
	{
		const double s = math::randu();
		const double t = 2 * M_PI * (i + s) / nn;
		xn[5 * i + 2] = t;
		xn[5 * i + 4] = sin(t);
	}
	for(uint32_t i = 0; i < nr; i++)
	{
		const double s = math::randu();
		const double t = 2 * M_PI * (i + s) / nr;
		xr[6 * i + 0] = t;
		xr[6 * i + 3] = sin(t);
	}
	//validator
	math::validation::Validator validator;
	validator.create_item();
	validator.item(0)->load_numeric(xn, nn, 5, 2, 4);
	validator.item(0)->load_reference(xr, nr, 6, 0, 3);
	//validate
	validator.validate();
}