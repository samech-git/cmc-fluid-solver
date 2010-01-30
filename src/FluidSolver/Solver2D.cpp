#include "Solver2D.h"

namespace FluidSolver
{
	void Solver2D::GetLayer(Vec2D *v, double *T, int outdimx, int outdimy)
	{
		if (outdimx == 0) outdimx = dimx;
		if (outdimy == 0) outdimy = dimy;

		for (int i = 0; i < outdimx; i++)
			for (int j = 0; j < outdimy; j++)
			{
				int x = (i * dimx / outdimx);
				int y = (j * dimy / outdimy);
				v[i * outdimy + j].x = next->U(x, y); 
				v[i * outdimy + j].y = next->V(x, y);
				T[i * outdimy + j] = next->T(x, y);
			}
	}

	void Solver2D::SetLayer(Vec2D *v, double *T)
	{
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				{
					cur->U(i, j) = v[i * dimy + j].x;
					cur->V(i, j) = v[i * dimy + j].y;
					cur->T(i, j) = T[i * dimy + j];
				}		
	}

	void Solver2D::UpdateBoundaries()
	{
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				switch(grid->GetType(i, j))
				{
				case BOUND:
				case VALVE:
					cur->U(i, j) = grid->GetData(i, j).vel.x;
					cur->V(i, j) = grid->GetData(i, j).vel.y;
					cur->T(i, j) = grid->GetData(i, j).T;
					break;
				}
		cur->CopyAllto(grid, next, BOUND);
		cur->CopyAllto(grid, next, VALVE);
	}

	void Solver2D::ReturnBoundaries()
	{
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
			{
				CondData2D d = grid->GetData(i, j);
				grid->SetFieldData(i, j, CondData2D(d.type, d.cell, Vec2D(cur->U(i, j), cur->V(i, j)), 0));
			}
	}
}