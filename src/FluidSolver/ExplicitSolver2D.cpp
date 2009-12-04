#include "ExplicitSolver2D.h"

namespace FluidSolver
{
	void ExplicitSolver2D::Init(Grid2D &_grid, FluidParams &_params)
	{
		FreeMemory();
	
		dimx = _grid.dimx;
		dimy = _grid.dimy;

		grid = new Grid2D(_grid);
		params = _params;

		cur = new TimeLayer2D(grid->dimx, grid->dimy, grid->dx, grid->dy);
		temp = new TimeLayer2D(grid->dimx, grid->dimy, grid->dx, grid->dy);
		next = new TimeLayer2D(grid->dimx, grid->dimy, grid->dx, grid->dy);

		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				switch(grid->GetType(i, j))
				{
				case IN:
				case BOUND:
				case VALVE:
				case OUT:
					cur->U(i, j) = grid->GetData(i, j).vel.x;
					cur->V(i, j) = grid->GetData(i, j).vel.y;
					cur->T(i, j) = grid->GetData(i, j).T;
					break;
				}

		cur->CopyAllto(grid, next);
		cur->CopyAllto(grid, temp);
	}

	void ExplicitSolver2D::SolveU(double dt, int num_local, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next)
	{
		TimeLayer2D *next_loc = new TimeLayer2D(grid->dimx, grid->dimy, grid->dx, grid->dy);

		for (int it = 0; it < num_local; it++)
		{
			// eval new layer
			for (int i = 0; i < dimx; i++)
				for (int j = 0; j < dimy; j++)
					if (grid->GetType(i, j) == IN)
					{
						next_loc->U(i, j) = cur->U(i, j) + dt * ( 
							- temp->U(i, j) * temp->Ux(i, j) 
							- temp->V(i, j) * temp->Uy(i, j) 
							- temp->Tx(i, j)
							+ (1 / params.Re) * (temp->Uxx(i, j) + temp->Uyy(i, j)) );
					}
			next_loc->MergeUto(grid, next, IN);
			cur->CopyUto(grid, next, BOUND);
			cur->CopyUto(grid, next, VALVE);
		}

		delete next_loc;
	}

	void ExplicitSolver2D::SolveV(double dt, int num_local, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next)
	{
		TimeLayer2D *next_loc = new TimeLayer2D(grid->dimx, grid->dimy, grid->dx, grid->dy);

		for (int it = 0; it < num_local; it++)
		{
			// eval new layer
			for (int i = 0; i < dimx; i++)
				for (int j = 0; j < dimy; j++)
					if (grid->GetType(i, j) == IN)
					{
						next_loc->V(i, j) = cur->V(i, j) + dt * ( 
							- temp->U(i, j) * temp->Vx(i, j) 
							- temp->V(i, j) * temp->Vy(i, j) 
							- temp->Ty(i, j)
							+ (1 / params.Re) * (temp->Vxx(i, j) + temp->Vyy(i, j)) );
					}
			next_loc->MergeVto(grid, next, IN);
			cur->CopyVto(grid, next, BOUND);
			cur->CopyVto(grid, next, VALVE);
		}

		delete next_loc;
	}

	void ExplicitSolver2D::SolveT(double dt, int num_local, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next)
	{
		TimeLayer2D *next_loc = new TimeLayer2D(grid->dimx, grid->dimy, grid->dx, grid->dy);

		for (int it = 0; it < num_local; it++)
		{
			// eval new layer
			for (int i = 0; i < dimx; i++)
				for (int j = 0; j < dimy; j++)
					if (grid->GetType(i, j) == IN)
					{
						next_loc->T(i, j) = cur->T(i, j) + dt * ( 
							- temp->U(i, j) * temp->Tx(i, j) 
							- temp->V(i, j) * temp->Ty(i, j) 
							+ (1 / (params.Re * params.Pr)) * (temp->Txx(i, j) + temp->Tyy(i, j)) 
							+ ((params.lambda - 1) / (params.lambda * params.Re)) * temp->DissFunc(i, j));
					}
			next_loc->MergeTto(grid, next, IN);
			cur->CopyTto(grid, next, BOUND);
			cur->CopyTto(grid, next, VALVE);
		}

		delete next_loc;
	}

	void ExplicitSolver2D::TimeStep(double dt, int num_global, int num_local)
	{
		// do global iterations
		int it;
		double err = next->EvalDivError(grid);

		for (it = 0; (it < num_global) || (err > ERR_THRESHOLD); it++)
		{
			SolveU(dt, num_local, cur, temp, next);
			SolveV(dt, num_local, cur, temp, next);
			SolveT(dt, num_local, cur, temp, next);

			err = next->EvalDivError(grid);

			next->MergeUto(grid, temp, IN);
			next->MergeVto(grid, temp, IN);
			next->MergeTto(grid, temp, IN);

			if (it > MAX_GLOBAL_ITERS) { printf("Exceeded max number of iterations (%i)\n", MAX_GLOBAL_ITERS); exit(1); }
		}

		// output number of global iterations
		printf("%i - ", it);

		next->CopyAllto(grid, cur);
	}

	void ExplicitSolver2D::GetResult(int nx, int ny, Vec2D *vel, double *T)
	{
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				int x = (i * dimx / nx);
				int y = (j * dimy / ny);
				vel[i * ny + j].x = next->U(x, y); 
				vel[i * ny + j].y = next->V(x, y);
				T[i * ny + j] = next->T(x, y);
			}
	}

	void ExplicitSolver2D::FreeMemory()
	{
		if (grid != NULL) delete grid;

		if (cur != NULL) delete cur;
		if (temp != NULL) delete temp;
		if (next != NULL) delete next;
	}

	ExplicitSolver2D::ExplicitSolver2D()
	{
		grid = NULL;
		cur = NULL;
		temp = NULL;
		next = NULL;
	}

	ExplicitSolver2D::~ExplicitSolver2D()
	{
		FreeMemory();
	}
}