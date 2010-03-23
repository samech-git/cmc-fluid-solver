#include "FluidSolver2D.h"

using namespace FluidSolver2D;

int main(int argc, char **argv)
{
	char inputPath[MAX_PATH];
	char outputPath[MAX_PATH];
	char configPath[MAX_PATH];

	FindFile(inputPath, argv[1]);
	FindFile(outputPath, argv[2], false);
	FindFile(configPath, argv[3]);
	
	Config::Config();
	Config::LoadFromFile(configPath);

	//--------------------------------------- Initializing ---------------------------------------
	Grid2D grid(Config::dx, Config::dy, Config::startT, Config::bc_noslip);
	if (grid.LoadFromFile(inputPath))
	{
		printf("dx,dy,dimx,dimy,bc_noslip\n");
		printf("%f,%f,%i,%i,%i\n", Config::dx, Config::dy, grid.dimx, grid.dimy, Config::bc_noslip);
	}
	grid.Prepare(0, 0);
	
	//FluidParams params(Re, Pr, lambda);
	FluidParams params(Config::viscosity, Config::density, Config::R_specific, Config::k, Config::cv);

	Solver2D *solver;
	switch (Config::solverID)
	{
		case Explicit: solver = new ExplicitSolver2D(); break;
		case ADI: solver = new AdiSolver2D(); break;
		case Stable: solver = new StableSolver2D(); break;
	}
	solver->Init(&grid, params);

	printf("Starting from the beginning\n");
	int startFrame = 0;
	FILE *resFile = NULL;
	fopen_s(&resFile, outputPath, "w");
	OutputResultHeader(resFile, &grid.bbox, Config::outdimx, Config::outdimy);
	
	Vec2D *resVel = new Vec2D[Config::outdimx * Config::outdimy];
	double *resT = new double[Config::outdimx * Config::outdimy];

	//------------------------------------------ Solving ------------------------------------------
	cpu_timer timer;
	timer.start();

	int frames = grid.GetFramesNum();
	double length = grid.GetCycleLenght();
	double dt = length / (frames * Config::calc_subframes);
	double finaltime = length * Config::cycles;

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
		solver->TimeStep(dt, Config::num_global, Config::num_local);
		solver->SetGridBoundaries();

		timer.stop();
		PrintTimeStepInfo(curentframe, i, t, finaltime, timer.elapsed_sec());

		if ((i % Config::out_subframes) == 0)
		{
			float dur = (float)dt * Config::out_subframes;
			if (dur > layer_time) dur = layer_time;

			solver->GetLayer(resVel, resT, Config::outdimx, Config::outdimy);
			OutputResult(resFile, resVel, resT, Config::outdimx, Config::outdimy, dur);
		}
	}
	printf("\n");

	delete solver;
	delete [] resVel;
	delete [] resT;

	fclose(resFile);
	return 0;
}