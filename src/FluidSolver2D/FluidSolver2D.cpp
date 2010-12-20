#include "FluidSolver2D.h"
#include <string.>

using namespace FluidSolver2D;

int main(int argc, char **argv)
{
	printf("%s\n", argv[0]);
	char inputPath[MAX_PATH];
	char outputPath[MAX_PATH];
	char curOutFile[MAX_PATH];
	char configPath[MAX_PATH];
	char fieldPath[MAX_PATH];

	if (argc == 2)
		LoadProject(argv[1], inputPath, fieldPath, outputPath, configPath);
	else
		if (argc > 4)
		{
			FindFile(inputPath, argv[1]);
			FindFile(fieldPath, argv[2]);
			FindFile(outputPath, argv[3], false);
			FindFile(configPath, argv[4]);
		}
		else
		{
			FindFile(inputPath, argv[1]);
			FindFile(outputPath, argv[2], false);
			FindFile(configPath, argv[3]);
			sprintf_s(fieldPath, MAX_STR_SIZE, "");
		}

	Config::Config();
	Config::LoadFromFile(configPath);

	//--------------------------------------- Initializing ---------------------------------------
	Grid2D grid(Config::dx, Config::dy, Config::startT, Config::bc_noslip, Config::bc_strength);
	if (grid.LoadFromFile(inputPath, fieldPath))
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

	Vec2D *resVel = new Vec2D[Config::outdimx * Config::outdimy];
	double *resT = new double[Config::outdimx * Config::outdimy];

	//------------------------------------------ Solving ------------------------------------------
	cpu_timer timer;
	timer.start();

	int frames = grid.GetFramesNum();
	double length = grid.GetCycleLenght();
	double dt = length / (frames * Config::calc_subframes);
	double finaltime = length * Config::cycles;

	OutputNetCDFHeader2D(outputPath, &grid.bbox, dt * Config::out_subframes, finaltime, Config::outdimx, Config::outdimy);

	printf("dt = %f\n", dt);
	int lastframe = -1;
	int currentcycle = 0;
	double t = dt;
	for (int i=0; t < finaltime; t+=dt, i++)
	{
		int currentframe = grid.GetFrame(t);
		float layer_time = grid.GetLayerTime(t);

		if (currentframe != lastframe)
		{
			if (currentframe == 0) //new cycle
			{
				currentcycle++;

				if( Config::outfmt == MultiVox )
				{
					if (currentcycle > 0)
					{
						char add[MAX_PATH];
						sprintf_s(add, MAX_STR_SIZE, "_%i", currentcycle);
						ExtendFileName(outputPath, curOutFile, add);
					}
					
					OutputResultHeader(curOutFile, &grid.bbox, Config::outdimx, Config::outdimy);
				}
			}

			if( Config::outfmt == MultiVox )
			{
				FILE *resFile = NULL;
				fopen_s(&resFile, curOutFile, "a");
				fprintf(resFile, "Frame %i\n", currentframe);
				fclose(resFile);
			}

			lastframe = currentframe;
			i = 0;
		}

		grid.Prepare(t);
		solver->UpdateBoundaries();
		solver->TimeStep((FTYPE)dt, Config::num_global, Config::num_local);
		solver->SetGridBoundaries();

		timer.stop();
		PrintTimeStepInfo(currentframe, i, t, finaltime, timer.elapsed_sec());

		if ((i % Config::out_subframes) == 0)
		{
			float dur = (float)dt * Config::out_subframes;
			if (dur > layer_time) dur = layer_time;

			solver->GetLayer(resVel, resT, Config::outdimx, Config::outdimy);
				
			if( Config::outfmt == MultiVox )
				OutputResult(curOutFile, resVel, resT, Config::outdimx, Config::outdimy, dur);
			else
				OutputNetCDF2D_U(outputPath, resVel, resT, Config::outdimx, Config::outdimy,  
								(i + Config::out_subframes >= Config::calc_subframes) && (currentframe == frames-1) && (currentcycle == Config::cycles));
		}
	}
	printf("\n");

	delete solver;
	delete [] resVel;
	delete [] resT;
	
	return 0;
}