#include "FluidSolver3D.h"

using namespace FluidSolver3D;
using namespace Common;

int main(int argc, char **argv)
{
	char configPath[MAX_STR_SIZE];

	Config::Config();
	Config::LoadFromFile(configPath);

	Grid3D grid(Config::dx, Config::dy, Config::dz);

	return 0;
}