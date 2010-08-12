#include "Solver3D.h"

namespace FluidSolver3D
{
	void Solver3D::GetLayer(Vec3D *v, double *T, int outdimx, int outdimy, int outdimz)
	{
		next->FilterToArrays(v, T, outdimx, outdimy, outdimz);
	}

	void Solver3D::UpdateBoundaries()
	{
		cur->CopyFromGrid(grid, NODE_BOUND);
		cur->CopyLayerTo(grid, next, NODE_BOUND);
	}

	void Solver3D::SetGridBoundaries()
	{
		cur->CopyToGrid(grid);
	}

	void Solver3D::ClearOutterCells()
	{
		next->Clear(grid, NODE_OUT, 0.0, 0.0, 0.0, (FTYPE)grid->baseT);
	}
}