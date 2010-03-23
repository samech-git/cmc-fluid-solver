#include "AdiSolver3D.h"

namespace FluidSolver3D
{
	AdiSolver3D::AdiSolver3D()
	{
		// TODO
	}

	AdiSolver3D::~AdiSolver3D()
	{
		// TODO
	}

	void AdiSolver3D::Init(Grid3D* _grid, FluidParams &_params)
	{
		grid = _grid;

		dimx = grid->dimx;
		dimy = grid->dimy;
		dimz = grid->dimz;

		params = _params;

		// TODO
	}

	void AdiSolver3D::TimeStep(double dt, int num_global, int num_local)
	{
		// TODO
	}
}