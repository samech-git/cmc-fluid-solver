#pragma once

#include "Grid2D.h"
#include "TimeLayer2D.h"

namespace FluidSolver2D
{
	class Solver2D
	{
	public:
		virtual void Init(Grid2D* grid, FluidParams &params) = 0;
		virtual void TimeStep(FTYPE dt, int num_global, int num_local) = 0;
		
		void GetLayer(Vec2D *v, double *T, int outdimx = 0, int outdimy = 0);
		void SetLayer(Vec2D *v, double *T);

		void UpdateBoundaries();
		void SetGridBoundaries();
		void ClearOutterCells();

		Grid2D *grid;

	protected:
		int dimx, dimy;
		FluidParams params;
		TimeLayer2D *cur, *next;

		double EvalDivError(TimeLayer2D *cur);
	};
}