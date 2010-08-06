#pragma once

#include "Grid3D.h"
#include "TimeLayer3D.h"

namespace FluidSolver3D
{
	class Solver3D
	{
	public:
		virtual void Init(Grid3D* grid, FluidParams &params) = 0;
		virtual void TimeStep(double dt, int num_global, int num_local) = 0;
		
		void GetLayer(Vec3D *v, double *T, int outdimx = 0, int outdimy = 0, int outdimz = 0);
		void SetLayer(Vec3D *v, double *T);

		void UpdateBoundaries();
		void SetGridBoundaries();
		void ClearOutterCells();

		Grid3D *grid;

		virtual ~Solver3D() {};

	protected:
		int dimx, dimy, dimz;
		FluidParams params;
		TimeLayer3D *cur, *next;

		double EvalDivError(TimeLayer3D *cur);
	};
}