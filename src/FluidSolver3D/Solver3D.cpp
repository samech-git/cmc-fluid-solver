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
					int ind = i * outdimy * outdimz + j * outdimz + k;
					v[ind].x = next->U->elem(x, y, z); 
					v[ind].y = next->V->elem(x, y, z);
					v[ind].z = next->W->elem(x, y, z);
					T[ind] = next->T->elem(x, y, z);
				}
	}

	void Solver3D::SetLayer(Vec3D *v, double *T)
	{
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				for (int k = 0; k < dimz; k++)
				{
					int ind = i * dimy * dimz + j * dimz + k;
					cur->U->elem(i, j, k) = v[ind].x;
					cur->V->elem(i, j, k) = v[ind].y;
					cur->W->elem(i, j, k) = v[ind].z;
					cur->T->elem(i, j, k) = T[ind];
				}		
	}

	void Solver3D::UpdateBoundaries()
	{
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				for (int k = 0; k < dimz; k++)
					if (grid->GetType(i, j, k) == NODE_BOUND)
					{
						Vec3D velocity = grid->GetVel(i, j, k);
						cur->U->elem(i, j, k) = velocity.x;
						cur->V->elem(i, j, k) = velocity.y;
						cur->W->elem(i, j, k) = velocity.z;
						cur->T->elem(i, j, k) = grid->GetT(i, j, k);
					}
		cur->CopyLayerTo(grid, next, NODE_BOUND);
	}

	void Solver3D::SetGridBoundaries()
	{
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				for (int k = 0; k < dimz; k++)
					grid->SetNodeVel(i, j, k, Vec3D(cur->U->elem(i, j, k), cur->V->elem(i, j, k), cur->W->elem(i, j , k)));
	}

	void Solver3D::ClearOutterCells()
	{
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				for (int k = 0; k < dimz; k++)
					if (grid->GetType(i, j, k) == NODE_OUT)
					{
						next->U->elem(i, j, k) = 0.0;
						next->V->elem(i, j, k) = 0.0;
						next->W->elem(i, j, k) = 0.0;
						next->T->elem(i, j, k) = grid->baseT;
					}
	}
}