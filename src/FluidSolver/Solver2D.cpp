#include "Solver2D.h"

namespace FluidSolver
{
	void Solver2D::GetResult(int outdimx, int outdimy, Vec2D *v, double *T)
	{
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
}