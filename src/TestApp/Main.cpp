#include "FluidSolver.h"

const double dx = 2.0;
const double dy = 2.0;

using namespace FluidSolver;

int main(int argc, char **argv)
{
	Test();

	Grid2D grid(dx, dy);
	if (grid.LoadFromFile("..\\..\\data\\test.txt") == OK)
		grid.TestPrint();

	return 0;
}