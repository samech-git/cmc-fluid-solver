#pragma once

#include "Grid2D.h"
#include "TimeLayer2D.h"

namespace FluidSolver
{
	class Solver2D
	{
	public:
		virtual void Init(Grid2D* grid, FluidParams &params) = 0;
		virtual void TimeStep(double dt, int num_global, int num_local) = 0;
		virtual void GetResult(int outdimx, int outdimy, Vec2D *v, double *T) = 0;

	protected:
		int dimx, dimy;

		double EvalDivError(TimeLayer2D *cur);
	};
}