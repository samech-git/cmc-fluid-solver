#include "FluidSolver.h"

const double dx = 0.25;
const double dy = 0.25;

const double dt = 0.1;

const double Re = 100.0;
const double Pr = 0.82;
const double lambda = 1.4;

const int num_global = 4;
const int num_local = 1;

const int nt = 1000;

const int outdimx = 50;
const int outdimy = 50;

using namespace FluidSolver;

int main(int argc, char **argv)
{
	Grid2D grid(dx, dy);
	if (grid.LoadFromFile("..\\..\\data\\test.txt") == OK)
	{
		printf("dx,dy,dimx,dimy,dt,Re,Pr,lambda\n");
		printf("%f,%f,%i,%i,%.3f,%f,%f,%f\n", dx, dy, grid.dimx, grid.dimy, dt, Re, Pr, lambda);
		//grid.TestPrint();
	}
	
	FluidParams params(Re, Pr, lambda);

	ExplicitSolver2D solver;
	solver.Init(grid, params);

	printf("\nglobal,local,err,time\n");
	for (int i = 0; i < nt; i++)
	{
 		solver.TimeStep(dt, num_global, num_local);
		printf("%.3f\n", i * dt);
	}

	Vec2D *vel = new Vec2D[outdimx * outdimy];
	double *T = new double[outdimx * outdimy];
	solver.GetResult(outdimx, outdimy, vel, T);

	TestPrintResult(outdimx, outdimy, vel, T);

	return 0;
}