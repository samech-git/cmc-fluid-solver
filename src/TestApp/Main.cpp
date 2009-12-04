#include "FluidSolver.h"

const double dx = 1.0;
const double dy = 1.0;

const double dt = 0.1;

const double Re = 10.0;
const double Pr = 0.82;
const double lambda = 1.4;

const int num_global = 2;
const int num_local = 2;

const int nt = 10000;

const int outdimx = 50;
const int outdimy = 50;

using namespace FluidSolver;

int main(int argc, char **argv)
{
	Grid2D grid(dx, dy);
	if (grid.LoadFromFile("..\\..\\data\\test.txt") == OK)
	{
		printf("%i %i\n", grid.dimx, grid.dimy);
		//grid.TestPrint();
	}
	
	FluidParams params(Re, Pr, lambda);

	ExplicitSolver2D solver;
	solver.Init(grid, params);

	for (int i = 0; i < nt; i++)
	{
		printf("%.3f\n", i * dt);
 		solver.TimeStep(dt, num_global, num_local);
	}

	Vec2D *vel = new Vec2D[outdimx * outdimy];
	double *T = new double[outdimx * outdimy];
	solver.GetResult(outdimx, outdimy, vel, T);

	TestPrintResult(outdimx, outdimy, vel, T);

	return 0;
}