#pragma once

#include <stdio.h>

namespace FluidSolver
{
	enum RetStatus { OK, ERR };

	static void Test()
	{
#ifdef _DEBUG
		printf("DEBUG\n");
#endif
		printf("Test!\n");
	}
}
