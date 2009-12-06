#include "FluidSolver.h"

const double dx = 0.5;
const double dy = 0.5;

const double dt = 0.5;

const double Re = 20.0;
const double Pr = 0.82;
const double lambda = 1.4;

const int num_global = 2;
const int num_local = 1;

//const int nt = 1000;
const int frames = 22;
const int subframes = 150;

const int outdimx = 50;
const int outdimy = 50;

using namespace FluidSolver;

int main(int argc, char **argv)
{
	//--------------------------------------- Initializing ---------------------------------------
	Grid2D grid(dx, dy);
	if (grid.LoadFromFile("..\\..\\data\\test.txt") == OK)
	{
		printf("dx,dy,dimx,dimy,dt,Re,Pr,lambda\n");
		printf("%f,%f,%i,%i,%.3f,%f,%f,%f\n", dx, dy, grid.dimx, grid.dimy, dt, Re, Pr, lambda);
		//grid.TestPrint();
	}
	grid.Prepare(0);
	
	FluidParams params(Re, Pr, lambda);

	ExplicitSolver2D solver;
	//solver.grid = &grid;
	solver.Init(&grid, params);

	//------------------------------------------ Solving ------------------------------------------
	Vec2D *vel = new Vec2D[outdimx * outdimy];
	double *T = new double[outdimx * outdimy];

	FILE *file = NULL;
	fopen_s(&file, "results.txt", "w");
	fprintf(file, "%.2f %.2f %.2f %.2f\n", grid.bbox.pMin.x, grid.bbox.pMin.y, grid.bbox.pMax.x, grid.bbox.pMax.y);

	float ddx = (grid.bbox.pMax.x - grid.bbox.pMin.x) / outdimx;
	float ddy = (grid.bbox.pMax.y - grid.bbox.pMin.y) / outdimy;
	fprintf(file, "%.2f %.2f %i %i\n", ddx, ddy, outdimx, outdimy);
	fprintf(file, "22\n");

	int percent = frames * subframes;
	int step = 0;

	for (int i = 0; i < frames; i++)
	{
		grid.Prepare(i);
		solver.UpdateBoundaries();
		
		fprintf(file, "0.035\n");
		for (int j = 0; j < subframes; j++)
		{
 			solver.TimeStep(dt, num_global, num_local);
			printf(" frame %i, substep %i (%i\%)\n", i, j, step / percent);
			step += 100;
		}
		solver.GetResult(outdimx, outdimy, vel, T);
		ShiferTestPrintResult(outdimx, outdimy, vel, T, file);
	}

	fclose(file);
	return 0;
}