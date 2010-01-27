#include "FluidSolver.h"
#include "Timer.h"

const double dx = 0.0014;
const double dy = 0.0014;

const double dt = 0.00001;

const double Re = 10.0;
const double Pr = 0.82;
const double lambda = 1.4;

const int num_global = 2;
const int num_local = 1;

const int frames = 25;
const int subframes = 100;

const int outdimx = 50;
const int outdimy = 50;

const float timeValue = 0.035f;

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
	Vec2D *lastVel = new Vec2D[grid.dimx * grid.dimy];
	double *lastT = new double[grid.dimx * grid.dimy];
	int startFrame = LoadLastLayer(lastPath, lastVel, lastT, grid.dimx, grid.dimy, frames);

	FILE *resFile = NULL;
	if (startFrame == 0)
	{
		printf("Starting from the beginning\n");
		fopen_s(&resFile, resPath, "w");
		OutputResultHeader(resFile, &grid.bbox, outdimx, outdimy, frames);
	}
	else if (startFrame == frames)
	{
		printf("All done!\n");
		return 0;
	}
	else
	{
		printf("Starting from frame %i\n", startFrame);
		fopen_s(&resFile, resPath, "a");
		solver->SetLayer(lastVel, lastT);
	}

	Vec2D *resVel = new Vec2D[outdimx * outdimy];
	double *resT = new double[outdimx * outdimy];

	//------------------------------------------ Solving ------------------------------------------
	for (int i = startFrame; i < frames; i++)
	{	
		for (int j = 0; j < subframes; j++)
		{
			cpu_timer timer;
			timer.start();

			grid.Prepare(i, (double)j / subframes);
			solver->UpdateBoundaries();
 			solver->TimeStep(dt, num_global, num_local);

			timer.stop();
			PrintTimeStepInfo(i, j, frames, subframes, timer.elapsed_sec());
		}
		
		solver->GetLayer(resVel, resT, outdimx, outdimy);
		OutputResult(resFile, resVel, resT, outdimx, outdimy, timeValue);

		solver->GetLayer(lastVel, lastT);
		SaveLastLayer(lastPath, i+1, lastVel, lastT, grid.dimx, grid.dimy); 
	}

	delete solver;
	delete [] resVel;
	delete [] resT;
	delete [] lastVel;
	delete [] lastT;

	fclose(resFile);
	return 0;
}