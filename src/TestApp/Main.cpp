#include "FluidSolver.h"
#include "Timer.h"

const double dx = 0.5;
const double dy = 0.5;

const double dt = 0.5;

const double Re = 15.0;
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
	if (grid.LoadFromFile("..\\..\\data\\test_heart.txt") == OK)
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
			cpu_timer timer;
			timer.start();

			grid.Prepare(i, (double)j / subframes);
			solver->UpdateBoundaries();
 			solver->TimeStep(dt, num_global, num_local);

			timer.stop();
			float time_left_sec = ((frames-i-1) * subframes + subframes-j-1) * timer.elapsed_sec();
			int time_h = ((int)time_left_sec) / 3600;
			int time_m = (((int)time_left_sec) / 60) % 60;
			int time_s = ((int)time_left_sec) % 60;
			printf(" frame %i\tsubstep %i\t%i%%\t(%i h %i m %i s left)\n", i, j, step / percent, time_h, time_m, time_s);
			step += 100;
		}
		solver->GetResult(outdimx, outdimy, vel, T);
		ShiferTestPrintResult(outdimx, outdimy, vel, T, file);
	}

	delete solver;

	fclose(file);
	return 0;
}