#include "FluidSolver.h"

const double dx = 0.5;
const double dy = 0.5;

const double dt = 1.0;

const double Re = 10.0;
const double Pr = 0.82;
const double lambda = 1.4;

const int num_global = 2;
const int num_local = 1;

const int frames = 22;
const int subframes = 150;

const int outdimx = 50;
const int outdimy = 50;

enum solvers { Explicit, ADI };
const int solverID = ADI;		

using namespace FluidSolver;

int main(int argc, char **argv)
{
	//--------------------------------------- Initializing ---------------------------------------
	Grid2D grid(dx, dy);
	if (grid.LoadFromFile("..\\..\\data\\test3.txt") == OK)
	{
		printf("dx,dy,dimx,dimy,dt,Re,Pr,lambda\n");
		printf("%f,%f,%i,%i,%.3f,%f,%f,%f\n", dx, dy, grid.dimx, grid.dimy, dt, Re, Pr, lambda);
	}
	grid.Prepare(0, 0);
	grid.TestPrint();
	
	FluidParams params(Re, Pr, lambda);

	Solver2D *solver;
	switch (solverID)
	{
		case Explicit: solver = new ExplicitSolver2D(); break;
		case ADI: solver = new AdiSolver2D(); break;
	}
	solver->Init(&grid, params);

	//------------------------------------------ Solving ------------------------------------------
	Vec2D *vel = new Vec2D[outdimx * outdimy];
	double *T = new double[outdimx * outdimy];

	FILE *file = NULL;
	fopen_s(&file, "results.txt", "w");
	fprintf(file, "%.2f %.2f %.2f %.2f\n", grid.bbox.pMin.x, grid.bbox.pMin.y, grid.bbox.pMax.x, grid.bbox.pMax.y);

	float ddx = (float)(grid.bbox.pMax.x - grid.bbox.pMin.x) / outdimx;
	float ddy = (float)(grid.bbox.pMax.y - grid.bbox.pMin.y) / outdimy;
	fprintf(file, "%.2f %.2f %i %i\n", ddx, ddy, outdimx, outdimy);
	fprintf(file, "%i\n", frames);

	int percent = frames * subframes;
	int step = 0;

	for (int i = 0; i < frames; i++)
	{	
		fprintf(file, "0.035\n");
		for (int j = 0; j < subframes; j++)
		{
			grid.Prepare(i, (double)j / subframes);
			solver->UpdateBoundaries();

 			solver->TimeStep(dt, num_global, num_local);
			printf(" frame %i\tsubstep %i\t%i%%\n", i, j, step / percent);
			step += 100;
		}
		solver->GetResult(outdimx, outdimy, vel, T);
		ShiferTestPrintResult(outdimx, outdimy, vel, T, file);
	}

	delete solver;

	fclose(file);
	return 0;
}