#include "FluidSolver.h"
#include "Timer.h"

const double dx = 0.25;
const double dy = 0.25;

const double dt = 0.15;

const double Re = 30.0;
const double Pr = 0.82;
const double lambda = 1.4;

const int num_global = 2;
const int num_local = 1;

const int frames = 25;
const int subframes = 150;

const int outdimx = 50;
const int outdimy = 50;

enum solvers { Explicit, ADI };
const int solverID = ADI;		

using namespace FluidSolver;

int main(int argc, char **argv)
{
	char dataPath[MAX_PATH];
	char resPath[MAX_PATH];
	char lastPath[MAX_PATH];
	
	sprintf_s(dataPath, "..\\..\\data\\%s_ns.txt", argv[1]);
	sprintf_s(resPath, "..\\..\\data\\%s_res.txt", argv[1]);
	sprintf_s(lastPath, "..\\..\\data\\%s_layer.txt", argv[1]);

	//--------------------------------------- Initializing ---------------------------------------
	Grid2D grid(dx, dy);
	if (grid.LoadFromFile(dataPath) == OK)
	{
		printf("dx,dy,dimx,dimy,dt,Re,Pr,lambda\n");
		printf("%f,%f,%i,%i,%.3f,%f,%f,%f\n", dx, dy, grid.dimx, grid.dimy, dt, Re, Pr, lambda);
	}
	grid.Prepare(0, 0);
	//grid.TestPrint();
	
	FluidParams params(Re, Pr, lambda);

	Solver2D *solver;
	switch (solverID)
	{
		case Explicit: solver = new ExplicitSolver2D(); break;
		case ADI: solver = new AdiSolver2D(); break;
	}
	solver->Init(&grid, params);

	// loading last-layer if there is one
	Vec2D *lastVel = NULL;
	double *lastT = NULL;
	int curFrame = LoadLastLayer(lastPath, &lastVel, &lastT, grid.dimx, grid.dimy);
	
	// result header
	Vec2D *resVel = new Vec2D[outdimx * outdimy];
	double *resT = new double[outdimx * outdimy];

	FILE *resFile = NULL;
	if (curFrame == 0)
	{
		fopen_s(&resFile, resPath, "w");
		fprintf(resFile, "%.2f %.2f %.2f %.2f\n", grid.bbox.pMin.x, grid.bbox.pMin.y, grid.bbox.pMax.x, grid.bbox.pMax.y);

		float ddx = (float)(grid.bbox.pMax.x - grid.bbox.pMin.x) / outdimx;
		float ddy = (float)(grid.bbox.pMax.y - grid.bbox.pMin.y) / outdimy;
		fprintf(resFile, "%.2f %.2f %i %i\n", ddx, ddy, outdimx, outdimy);
		fprintf(resFile, "%i\n", frames);
	}
	else
	{
		fopen_s(&resFile, resPath, "a");
	}

	//------------------------------------------ Solving ------------------------------------------
	int percent = frames * subframes;
	int step = 0;

	for (int i = 0; i < frames; i++)
	{	
		fprintf(resFile, "0.035\n");
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
		
		solver->GetLayer(resVel, resT, outdimx, outdimy);
		ShiferTestPrintResult(resFile, resVel, resT, outdimx, outdimy);

		solver->GetLayer(lastVel, lastT);
		
		FILE *lastFile = NULL;
		fopen_s(&lastFile, lastPath, "w");
		fprintf(lastFile, "%i %i\n", grid.dimx, grid.dimy);
		PrintLast(lastFile, lastVel, lastT, grid.dimx, grid.dimy);
		fclose(lastFile);
	}

	delete solver;
	delete [] resVel;
	delete [] resT;
	delete [] lastVel;
	delete [] lastT;

	fclose(resFile);

	return 0;
}