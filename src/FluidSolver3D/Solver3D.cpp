#include "Solver3D.h"

namespace FluidSolver3D
{
	void Solver3D::GetLayer(Vec3D *v, double *T, int outdimx, int outdimy, int outdimz)
	{
		if (outdimx == 0) outdimx = dimx;
		if (outdimy == 0) outdimy = dimy;
		if (outdimz == 0) outdimz = dimz;

		for (int i = 0; i < outdimx; i++)
			for (int j = 0; j < outdimy; j++)
				for (int k = 0; k < outdimz; k++)
				{	
					int x = (i * dimx / outdimx);
					int y = (j * dimy / outdimy);
					int z = (k * dimz / outdimz);
					v[i * outdimy * outdimz + j * outdimz + k].x = next->U->elem(x, y, z); 
					v[i * outdimy * outdimz + j * outdimz + k].y = next->V->elem(x, y, z);
					v[i * outdimy * outdimz + j * outdimz + k].z = next->W->elem(x, y, z);
					T[i * outdimy * outdimz + j * outdimz + k] = next->T->elem(x, y, z);
				}
	}

	void Solver3D::SetLayer(Vec3D *v, double *T)
	{
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				for (int k = 0; k < dimz; k++)
				{
					cur->U->elem(i, j, k) = v[i * dimy * dimz + j * dimz + k].x;
					cur->V->elem(i, j, k) = v[i * dimy * dimz + j * dimz + k].y;
					cur->W->elem(i, j, k) = v[i * dimy * dimz + j * dimz + k].z;
					cur->T->elem(i, j, k) = T[i * dimy * dimz + j * dimz + k];
				}		
	}

	void Solver3D::UpdateBoundaries()
	{
	}

	void Solver3D::SetGridBoundaries()
	{
	}

	void Solver3D::ClearOutterCells()
	{
	}
}