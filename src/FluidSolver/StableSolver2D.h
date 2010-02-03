#pragma once

#include "Solver2D.h"
#include "ScalarField2D.h"
#include <omp.h>

#define DIV_ERR_THRESHOLD		0.1
#define POISSON_ERR_THRESHOLD	1e-2
#define MAX_GLOBAL_ITERS		100

namespace FluidSolver
{
	class StableSolver2D : public Solver2D
	{
	public:
		void Init(Grid2D* _grid, FluidParams &_params);
		void TimeStep(double dt, int num_global, int num_local);

		StableSolver2D();
		~StableSolver2D();
	
	private:
		TimeLayer2D *temp, *next_w;
		ScalarField2D *q[2];
		ScalarField2D *div;

		void SolveU(double dt, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next);
		void SolveV(double dt, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next);
		void Project(double dt, TimeLayer2D *w, TimeLayer2D *proj);

		void FreeMemory();
	};
}