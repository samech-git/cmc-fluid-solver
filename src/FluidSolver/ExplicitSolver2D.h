#pragma once

#include "Solver2D.h"

#define ERR_THRESHOLD		0.1
#define MAX_GLOBAL_ITERS	100

namespace FluidSolver
{
	class ExplicitSolver2D : Solver2D
	{
	public:
		void Init(Grid2D* _grid, FluidParams &_params);
		void TimeStep(double dt, int num_global, int num_local);
		void GetResult(int outdimx, int outdimy, Vec2D *vel, double *T);

		ExplicitSolver2D();
		~ExplicitSolver2D();

		Grid2D *grid;
		void UpdateBoundaries();

	protected:
		//Grid2D *grid;
		TimeLayer2D *cur, *temp, *next;
		FluidParams params;
	
	private:
		void SolveU(double dt, int num_local, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next);
		void SolveV(double dt, int num_local, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next);
		void SolveT(double dt, int num_local, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next);

		void FreeMemory();
	};
}