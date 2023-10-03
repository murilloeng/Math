//maker
#include "../Maker/inc/Maker.hpp"

int main(int argc, char** argv)
{
	//setup
	Maker maker;
	maker.m_out = "Math";
	maker.setup(argc, argv);
	//build
	if(!maker.m_clean)
	{
		maker.build_src();
		maker.build_lib();
	}
	if(maker.m_clean)
	{
		maker.build_clean();
	}
	//return
	return EXIT_SUCCESS;
}