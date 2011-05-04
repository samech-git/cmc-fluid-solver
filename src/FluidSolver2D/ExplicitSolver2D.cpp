/*
 *  Copyright 2010-2011 Nikolai Sakharnykh
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include "ExplicitSolver2D.h"

namespace FluidSolver2D
{
	void ExplicitSolver2D::Init(Grid2D* _grid, FluidParams &_params)
	{
		grid = _grid;

		dimx = grid->dimx;
		dimy = grid->dimy;

		params = _params;

		cur = new TimeLayer2D(grid->dimx, grid->dimy, (FTYPE)grid->dx, (FTYPE)grid->dy);
		next = new TimeLayer2D(grid->dimx, grid->dimy, (FTYPE)grid->dx, (FTYPE)grid->dy);
		
		temp = new TimeLayer2D(grid->dimx, grid->dimy, (FTYPE)grid->dx, (FTYPE)grid->dy);
		next_local = new TimeLayer2D(grid->dimx, grid->dimy, (FTYPE)grid->dx, (FTYPE)grid->dy);

		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				switch(grid->GetType(i, j))
				{
				case NODE_IN:
				case NODE_BOUND:
				case NODE_VALVE:
				case NODE_OUT:
					cur->U(i, j) = grid->GetData(i, j).vel.x;
					cur->V(i, j) = grid->GetData(i, j).vel.y;
					cur->T(i, j) = grid->GetData(i, j).T;
					break;
				}

		cur->CopyAllto(grid, next);
		cur->CopyAllto(grid, temp);
	}

	void ExplicitSolver2D::SolveU(FTYPE dt, int num_local, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next)
	{
		for (int it = 0; it < num_local; it++)
		{
			// eval new layer
			for (int i = 0; i < dimx; i++)
				for (int j = 0; j < dimy; j++)
					if (grid->GetType(i, j) == NODE_IN)
					{
						next_local->U(i, j) = cur->U(i, j) + dt * ( 
							- temp->U(i, j) * temp->Ux(i, j) 
							- temp->V(i, j) * temp->Uy(i, j) 
							- params.v_T * temp->Tx(i, j)
							+ params.v_vis * (temp->Uxx(i, j) + temp->Uyy(i, j)) );
					}
			next_local->CopyUto(grid, next, NODE_IN);
		}
	}

	void ExplicitSolver2D::SolveV(FTYPE dt, int num_local, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next)
	{
		for (int it = 0; it < num_local; it++)
		{
			// eval new layer
			for (int i = 0; i < dimx; i++)
				for (int j = 0; j < dimy; j++)
					if (grid->GetType(i, j) == NODE_IN)
					{
						next_local->V(i, j) = cur->V(i, j) + dt * ( 
							- temp->U(i, j) * temp->Vx(i, j) 
							- temp->V(i, j) * temp->Vy(i, j) 
							- params.v_T * temp->Ty(i, j)
							+ params.v_vis * (temp->Vxx(i, j) + temp->Vyy(i, j)) );
					}
			next_local->CopyVto(grid, next, NODE_IN);
		}
	}

	void ExplicitSolver2D::SolveT(FTYPE dt, int num_local, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next)
	{
		for (int it = 0; it < num_local; it++)
		{
			// eval new layer
			for (int i = 0; i < dimx; i++)
				for (int j = 0; j < dimy; j++)
					if (grid->GetType(i, j) == NODE_IN)
					{
						next_local->T(i, j) = cur->T(i, j) + dt * ( 
							- temp->U(i, j) * temp->Tx(i, j) 
							- temp->V(i, j) * temp->Ty(i, j) 
							+ params.t_vis * (temp->Txx(i, j) + temp->Tyy(i, j))
							+ params.t_phi * temp->DissFunc(i, j));
					}
			next_local->CopyTto(grid, next, NODE_IN);
		}
	}

	void ExplicitSolver2D::TimeStep(FTYPE dt, int num_global, int num_local)
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
			next->MergeAllto(grid, temp, NODE_IN);
			
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

		ClearOutterCells();

		// output number of global iterations & error
		printf("\rerr = %.4f,", err);

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