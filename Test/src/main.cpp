//std
#include <cstdio>
#include <stdexcept>

//Test
#include "Math/Test/inc/fem.hpp"
#include "Math/Test/inc/groups.hpp"
#include "Math/Test/inc/solvers.hpp"
#include "Math/Test/inc/geometry.hpp"
#include "Math/Test/inc/rotations.hpp"
#include "Math/Test/inc/miscellaneous.hpp"

#include "Math/inc/solvers/Harmonic.hpp"

void test(void)
{
	//open
	uint32_t rows = 0, cols = 0;
	FILE* file = fopen("Test/data/solvers/harmonic/duffing/numeric.dat", "r");
	//rows
	while(!feof(file))
	{
		if(fgetc(file) == '\n') rows++;
	}
	//cols
	rewind(file);
	char line[2048];
	fgets(line, sizeof(line), file);
	for(char c : line)
	{
		if(c == ' ') cols++;
		else if(c == '\n') break;
	}
	printf("rows: %d cols: %d\n", rows, cols);
	//close
	fclose(file);
}

int main(void)
{
	try
	{
		test();
		// tests::solvers::harmonic::duffing();
	}
	catch(const std::exception& exception)
	{
		printf("%s\n", exception.what());
	}
	//return
	return EXIT_SUCCESS;
}