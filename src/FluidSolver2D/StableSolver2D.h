#pragma once

#include "Solver2D.h"
#include "ScalarField2D.h"
#include <omp.h>

#define DIV_ERR_THRESHOLD		0.1
#define POISSON_ERR_THRESHOLD	1e-2
#define MAX_GLOBAL_ITERS		100

namespace FluidSolver2D
{
	class StableSolver2D : public Solver2D
	{
	public:
		void Init(Grid2D* _grid, FluidParams &_params);
		void TimeStep(FTYPE dt, int num_global, int num_local);

		StableSolver2D();
		~StableSolver2D();
	
	private:
		TimeLayer2D *temp, *next_w;
		ScalarField2D *q[2];
		ScalarField2D *div;

		int numInner, numBound;
		int *innerIndices;
		int *boundIndices;

		void PrepareIndices();
		void SolveU(FTYPE dt, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next);
		void SolveV(FTYPE dt, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next);
		void Project(FTYPE dt, TimeLayer2D *w, TimeLayer2D *proj);

		void FreeMemory();
	};
}