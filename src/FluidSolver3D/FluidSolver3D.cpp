#include "FluidSolver3D.h"

using namespace FluidSolver3D;
using namespace Common;

int main(int argc, char **argv)
{
	char inputPath[MAX_STR_SIZE];
	char outputPath[MAX_STR_SIZE];
	char configPath[MAX_STR_SIZE];

	FindFile(inputPath, argv[1]);
	FindFile(outputPath, argv[2], false);
	FindFile(configPath, argv[3]);

	Config::Config();
	Config::LoadFromFile(configPath);

	//--------------------------------------- Initializing ---------------------------------------
	Grid3D grid(Config::dx, Config::dy, Config::dz, Config::depth);
	if (grid.LoadFromFile(inputPath))
	{
		printf("dx,dy,dz,dimx,dimy,dimz,bc_noslip\n");
		printf("%f,%f,%f,%i,%i,%i,%i\n", Config::dx, Config::dy, Config::dz, grid.dimx, grid.dimy, grid.dimz, Config::bc_noslip);
	}
	grid.Prepare(0.0);

	FluidParams params(Config::viscosity, Config::density, Config::R_specific, Config::k, Config::cv);

	Solver3D *solver;
	switch (Config::solverID)
	{
		case Explicit: printf("Explicit solver is not implemented yet!\n"); break;
		case ADI: solver = new AdiSolver3D(); break;
		case Stable: printf("Stable solver is not implemented yet!\n"); break;
	}
	solver->Init(&grid, params);


	return 0;
}