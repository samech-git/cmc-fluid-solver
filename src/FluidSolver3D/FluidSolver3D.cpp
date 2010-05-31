#include "FluidSolver3D.h"

using namespace FluidSolver3D;
using namespace Common;

void WriteNewFrame(vector<string> &outputPath, int currentframe)
{
	for (int z = 0; z < Config::outdimz; z++)
	{
		FILE *resFile = NULL;
		fopen_s(&resFile, outputPath[z].c_str(), "a");
		fprintf(resFile, "Frame %i\n", currentframe);
		fclose(resFile);
	}
}

int main(int argc, char **argv)
{
	char inputPath[MAX_STR_SIZE];
	char configPath[MAX_STR_SIZE];
	vector<string> outputPath;

	FindFile(inputPath, argv[1]);
	FindFile(configPath, argv[3]);

	Config::Config();
	Config::LoadFromFile(configPath);

	//--------------------------------------- Initializing ---------------------------------------
	Grid3D grid(Config::dx, Config::dy, Config::dz, Config::depth, Config::startT);
	if (grid.LoadFromFile(inputPath))
	{
		printf("dx,dy,dz,dimx,dimy,dimz,bc_noslip\n");
		printf("%f,%f,%f,%i,%i,%i,%i\n", Config::dx, Config::dy, Config::dz, grid.dimx, grid.dimy, grid.dimz, Config::bc_noslip);
	}
	grid.Prepare(0.0);

	grid.TestPrint("grid.txt");

	FluidParams *params;
	if (Config::useNormalizedParams) params = new FluidParams(Config::Re, Config::Pr, Config::lambda);
		else params = new FluidParams(Config::viscosity, Config::density, Config::R_specific, Config::k, Config::cv);

	Solver3D *solver;
	switch (Config::solverID)
	{
		case Explicit: printf("Explicit solver is not implemented yet!\n"); break;
		case ADI: solver = new AdiSolver3D(); break;
		case Stable: printf("Stable solver is not implemented yet!\n"); break;
	}
	solver->Init(&grid, *params);

	printf("Starting from the beginning\n");
	int startFrame = 0;
	
	for (int z = 0; z < Config::outdimz; z++)
	{
		char buf[MAX_STR_SIZE];
		sprintf_s(buf, MAX_STR_SIZE, "%s_%i.txt", argv[2], z);
		outputPath.push_back(buf);
		OutputResultHeader(outputPath[z].c_str(), &grid.GetGrid2D()->bbox, Config::outdimx, Config::outdimy);
	}

	// allocate result arrays
	Vec3D *resVel = new Vec3D[Config::outdimx * Config::outdimy * Config::outdimz];
	double *resT = new double[Config::outdimx * Config::outdimy * Config::outdimz];

	//------------------------------------------ Solving ------------------------------------------
	cpu_timer timer;
	timer.start();

	int frames = grid.GetGrid2D()->GetFramesNum();
	double length = grid.GetGrid2D()->GetCycleLenght();
	double dt = length / (frames * Config::calc_subframes);
	double finaltime = length * Config::cycles;

	printf("dt = %f\n", dt);
	int lastframe = -1;
	double t = dt;
	for (int i=0; t < finaltime; t+=dt, i++)
	{
		int currentframe = grid.GetGrid2D()->GetFrame(t);
		float layer_time = grid.GetGrid2D()->GetLayerTime(t);

		if (currentframe != lastframe)
		{
			WriteNewFrame(outputPath, currentframe);
			lastframe = currentframe;
			i = 0;
		}

		grid.Prepare(t);
		solver->UpdateBoundaries();
		solver->TimeStep(dt, Config::num_global, Config::num_local);
		solver->SetGridBoundaries();

		timer.stop();
		PrintTimeStepInfo(currentframe, i, t, finaltime, timer.elapsed_sec());

		if ((i % Config::out_subframes) == 0)
		{
			float dur = (float)dt * Config::out_subframes;
			if (dur > layer_time) dur = layer_time;

			solver->GetLayer(resVel, resT, Config::outdimx, Config::outdimy, Config::outdimz);
			for (int z = 0; z < Config::outdimz; z++)
				OutputSliceResult(outputPath[z].c_str(), z, resVel, resT, Config::outdimx, Config::outdimy, Config::outdimz, dur);
		}
	}
	printf("\n");

	delete solver;
	delete [] resVel;
	delete [] resT;

	return 0;
}