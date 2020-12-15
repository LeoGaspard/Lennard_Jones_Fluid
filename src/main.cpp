#include <iostream>
#include <exception>

#include "CError.h"
#include "CDynamic.h"
#include "Constants.h"

int main(int argc, const char * argv[])
{
	try
	{
		CDynamic *dyn = new CDynamic();
		dyn->Setup(argc, argv);
		dyn->Run();
	}
	catch(std::exception const& e)
	{
		std::cerr << "FATAL ERROR : " << e.what() << std::endl;
		exit(1);
	}
	return 0;
}
