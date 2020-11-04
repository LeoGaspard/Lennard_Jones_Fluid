#include <iostream>

#include "CDynamic.h"
#include "Constants.h"

int main(int argc, const char * argv[])
{
	CDynamic *dyn = new CDynamic();
	dyn->Setup(argc, argv);

	return 0;
}
