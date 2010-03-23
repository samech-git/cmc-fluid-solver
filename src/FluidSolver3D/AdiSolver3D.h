#pragma once

#include "Solver3D.h"

#define ERR_THRESHOLD		0.1

namespace FluidSolver3D
{
	class AdiSolver3D : public Solver3D
	{
	public:
		AdiSolver3D();
		~AdiSolver3D();

		void Init(Grid3D* _grid, FluidParams &_params);
		void TimeStep(double dt, int num_global, int num_local);
	};
}