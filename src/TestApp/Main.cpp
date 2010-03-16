#include "FluidSolver.h"
#include "Timer.h"
#define BC_NOSLIP 1

// grid size
const double dx = 0.0007;
const double dy = 0.0007;

// old params
const double Re = 50.0;
const double Pr = 0.82;
const double lambda = 1.4;

// new params
const double viscosity = 0.01;		// temporary high, 0.001002 for water at 20 C
const double density = 1000.0;		// water

// thermodynamic params
const double R_specific = 461.495;	// water,	287.058	for air (gas constant)
const double k = 0.6;				// water (thermal conductivity)
const double cv = 4200.0;			// water (specific heat capacity at constant volume)
const double startT = 300.0;		// in Kelvin 

// solver params
const int num_global = 2;
const int num_local = 1;

// animation params
const int cycles = 1;
const int subframes = 50;
const int out_subframes = 5;

// output grid
const int outdimx = 50;
const int outdimy = 50;

// solver type
enum solvers { Explicit, ADI, Stable };
const int solverID = Stable;		

using namespace FluidSolver;

int main(int argc, char **argv)
{
	char dataPath[MAX_PATH];
	char resPath[MAX_PATH];

	sprintf_s(dataPath, "..\\..\\data\\%s_ns.txt", argv[1]);
	sprintf_s(resPath, "..\\..\\data\\%s_res.txt", argv[1]);

	//--------------------------------------- Initializing ---------------------------------------
	Grid2D grid(dx, dy, startT);
	if (grid.LoadFromFile(dataPath) == OK)
	{
		printf("dx,dy,dimx,dimy,dt\n");
		printf("%f,%f,%i,%i\n", dx, dy, grid.dimx, grid.dimy);
	}
	grid.Prepare(0, 0);
	//grid.TestPrint();
	
	//FluidParams params(Re, Pr, lambda);
	FluidParams params(viscosity, density, R_specific, k, cv);

	Solver2D *solver;
	switch (solverID)
	{
		case Explicit: solver = new ExplicitSolver2D(); break;
		case ADI: solver = new AdiSolver2D(); break;
		case Stable: solver = new StableSolver2D(); break;
	}
	solver->Init(&grid, params);

	printf("Starting from the beginning\n");
	int startFrame = 0;
	FILE *resFile = NULL;
	fopen_s(&resFile, resPath, "w");
	OutputResultHeader(resFile, &grid.bbox, outdimx, outdimy);
	
	Vec2D *resVel = new Vec2D[outdimx * outdimy];
	double *resT = new double[outdimx * outdimy];

	//------------------------------------------ Solving ------------------------------------------
	cpu_timer timer;
	timer.start();

	int frames = grid.GetFramesNum();
	double length = grid.GetCycleLenght();
	double dt = length / (frames * subframes);
	double finaltime = length * cycles;

	printf("dt = %f\n", dt);
	int lastframe = -1;
	double t = dt;
	for (int i=0; t < finaltime; t+=dt, i++)
	{
		int curentframe = grid.GetFrame(t);
		float layer_time = grid.GetLayerTime(t);

		if (curentframe != lastframe)
		{
			fprintf(resFile, "Frame %i\n", curentframe);
			lastframe = curentframe;
			i = 0;
		}

		grid.Prepare(t);
		solver->UpdateBoundaries();
		solver->TimeStep(dt, num_global, num_local);
		solver->SetGridBoundaries();

		timer.stop();
		PrintTimeStepInfo(curentframe, i, t, finaltime, timer.elapsed_sec());

		if ((i % out_subframes) == 0)
		{
			float dur = dt * out_subframes;
			if (dur > layer_time) dur = layer_time;

			solver->GetLayer(resVel, resT, outdimx, outdimy);
			OutputResult(resFile, resVel, resT, outdimx, outdimy, dur);
		}
	}
	printf("\n");

	delete solver;
	delete [] resVel;
	delete [] resT;

	fclose(resFile);
	return 0;
}