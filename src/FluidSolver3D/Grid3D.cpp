#include "Grid3D.h"

namespace FluidSolver3D
{
	Grid3D::Grid3D(double _dx, double _dy, double _dz) : 
		dx(_dx), dy(_dy), dz(_dz), nodes(NULL)
	{
	}

	bool Grid3D::LoadFromFile(char *filename)
	{
		FILE *file = NULL;
		fopen_s(&file, filename, "r");
		if (!file) return false;



		fclose(file);
		return true;
	}
}