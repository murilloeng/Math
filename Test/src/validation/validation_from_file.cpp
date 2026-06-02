//std
#include <cmath>

//Math
#include "Math/Test/inc/validation.hpp"
#include "Math/inc/Miscellaneous/util.hpp"
#include "Math/inc/Validation/Validator.hpp"

void tests::validation::validation_from_file(void)
{
	//data
	srand(time(nullptr));
	const uint32_t nr = 100;
	const uint32_t nn = 100;
	FILE* file_numeric = fopen("numeric.txt", "w");
	FILE* file_reference = fopen("reference.txt", "w");
	for(uint32_t i = 0; i < nn; i++)
	{
		const double s = math::randu();
		const double t = 2 * M_PI * (i + s) / nn;
		fprintf(file_numeric, "%+.6e %+.6e \n", t, sin(t));
	}
	for(uint32_t i = 0; i < nr; i++)
	{
		const double s = math::randu();
		const double t = 2 * M_PI * (i + s) / nr;
		fprintf(file_reference, "%+.6e %+.6e \n", t, sin(t));
	}
	fclose(file_numeric);
	fclose(file_reference);
	//validator
	math::validation::Validator validator;
	validator.create_item();
	validator.item(0)->load_numeric("numeric.txt", 0, 1);
	validator.item(0)->load_reference("reference.txt", 0, 1);
	//validate
	validator.validate();
	//remove
	remove("numeric.txt");
	remove("reference.txt");
}