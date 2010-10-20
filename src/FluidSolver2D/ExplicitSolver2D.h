#pragma once

#include "Solver2D.h"

#define ERR_THRESHOLD		0.1
#define MAX_GLOBAL_ITERS	100

namespace FluidSolver2D
{
	class ExplicitSolver2D : public Solver2D
	{
	public:
		void Init(Grid2D* _grid, FluidParams &_params);
		void TimeStep(FTYPE dt, int num_global, int num_local);

		ExplicitSolver2D();
		~ExplicitSolver2D();
	
	private:
		TimeLayer2D *temp, *next_local;

		void SolveU(FTYPE dt, int num_local, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next);
		void SolveV(FTYPE dt, int num_local, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next);
		void SolveT(FTYPE dt, int num_local, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next);

		void FreeMemory();
	};
}