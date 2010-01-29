#include "ExplicitSolver2D.h"

namespace FluidSolver
{
	void ExplicitSolver2D::Init(Grid2D* _grid, FluidParams &_params)
	{
		grid = _grid;

		dimx = grid->dimx;
		dimy = grid->dimy;

		params = _params;

		cur = new TimeLayer2D(grid->dimx, grid->dimy, grid->dx, grid->dy);
		next = new TimeLayer2D(grid->dimx, grid->dimy, grid->dx, grid->dy);
		
		temp = new TimeLayer2D(grid->dimx, grid->dimy, grid->dx, grid->dy);
		next_local = new TimeLayer2D(grid->dimx, grid->dimy, grid->dx, grid->dy);

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
		for (int it = 0; it < num_local; it++)
		{
			// eval new layer
			for (int i = 0; i < dimx; i++)
				for (int j = 0; j < dimy; j++)
					if (grid->GetType(i, j) == IN)
					{
						next_local->U(i, j) = cur->U(i, j) + dt * ( 
							- temp->U(i, j) * temp->Ux(i, j) 
							- temp->V(i, j) * temp->Uy(i, j) 
							- temp->Tx(i, j)
							+ (1 / params.Re) * (temp->Uxx(i, j) + temp->Uyy(i, j)) );
					}
			next_local->CopyUto(grid, next, IN);
		}
	}

	void ExplicitSolver2D::SolveV(double dt, int num_local, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next)
	{
		for (int it = 0; it < num_local; it++)
		{
			// eval new layer
			for (int i = 0; i < dimx; i++)
				for (int j = 0; j < dimy; j++)
					if (grid->GetType(i, j) == IN)
					{
						next_local->V(i, j) = cur->V(i, j) + dt * ( 
							- temp->U(i, j) * temp->Vx(i, j) 
							- temp->V(i, j) * temp->Vy(i, j) 
							- temp->Ty(i, j)
							+ (1 / params.Re) * (temp->Vxx(i, j) + temp->Vyy(i, j)) );
					}
			next_local->CopyVto(grid, next, IN);
		}
	}

	void ExplicitSolver2D::SolveT(double dt, int num_local, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next)
	{
		for (int it = 0; it < num_local; it++)
		{
			// eval new layer
			for (int i = 0; i < dimx; i++)
				for (int j = 0; j < dimy; j++)
					if (grid->GetType(i, j) == IN)
					{
						next_local->T(i, j) = cur->T(i, j) + dt * ( 
							- temp->U(i, j) * temp->Tx(i, j) 
							- temp->V(i, j) * temp->Ty(i, j) 
							+ (1 / (params.Re * params.Pr)) * (temp->Txx(i, j) + temp->Tyy(i, j))
							+ ((params.lambda - 1) / (params.lambda * params.Re)) * temp->DissFunc(i, j));
					}
			next_local->CopyTto(grid, next, IN);
		}
	}

	void ExplicitSolver2D::TimeStep(double dt, int num_global, int num_local)
	{
		cur->CopyAllto(grid, temp);

		// do global iterations
		int it;
		double err = next->EvalDivError(grid);

		for (it = 0; (it < num_global) || (err > ERR_THRESHOLD); it++)
		{
			// solve equations
			SolveU(dt, num_local, cur, temp, next);
			SolveV(dt, num_local, cur, temp, next);
			SolveT(dt, num_local, cur, temp, next);

			err = next->EvalDivError(grid);
			
			// update non-linear parameters
			next->MergeAllto(grid, temp, IN);
			
			if (it > MAX_GLOBAL_ITERS) 
			{
				printf("\nExceeded max number of iterations (%i)\n", MAX_GLOBAL_ITERS); 
				exit(1); 
			}

			if (err > ERR_THRESHOLD * 10)
			{
				printf("\nError is too big!\n", err);
				exit(1);
			}
		}

		// clear outter cells
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				if (grid->GetType(i, j) == OUT)
				{
					next->U(i, j) = 0.0;
					next->V(i, j) = 0.0;
					next->T(i, j) = 1.0;
				}

		// output number of global iterations & error
		printf("\r%i,%i,%.4f,", it, num_local, err);

		// copy result to current layer
		next->CopyAllto(grid, cur);
	}

	void ExplicitSolver2D::FreeMemory()
	{
		if (cur != NULL) delete cur;
		if (temp != NULL) delete temp;
		if (next_local != NULL) delete next_local;
		if (next != NULL) delete next;
	}

	ExplicitSolver2D::ExplicitSolver2D()
	{
		grid = NULL;

		cur = NULL;
		temp = NULL;
		next_local = NULL;
		next = NULL;
	}

	ExplicitSolver2D::~ExplicitSolver2D()
	{
		FreeMemory();
	}
}