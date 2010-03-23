#include "Grid3D.h"

namespace FluidSolver3D
{
	Grid3D::Grid3D(double _dx, double _dy, double _dz, double _depth) : 
		dx(_dx), dy(_dy), dz(_dz), depth(_depth), nodes(NULL)
	{
		grid2D = new FluidSolver2D::Grid2D(dx, dy, 0, true);
	}

	Grid3D::~Grid3D()
	{
		if (nodes != NULL) delete [] nodes;
		if (grid2D != NULL) delete grid2D;
	}

	NodeType Grid3D::GetType(int i, int j, int k)
	{
		return nodes[i * dimy * dimz + j * dimz + k].type;
	}

	bool Grid3D::LoadFromFile(char *filename)
	{
		if (grid2D->LoadFromFile(filename))
		{
			dimx = grid2D->dimx;
			dimy = grid2D->dimy;
			dimz = (int)ceil(depth / dz) + 1;
			return true;
		}
		else
			return false;
	}

	void Grid3D::Prepare(double time)
	{
		grid2D->Prepare(time);
	}
}