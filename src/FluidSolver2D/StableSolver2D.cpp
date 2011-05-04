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

#include "StableSolver2D.h"

namespace FluidSolver2D
{
	void StableSolver2D::Init(Grid2D* _grid, FluidParams &_params)
	{
		grid = _grid;

		dimx = grid->dimx;
		dimy = grid->dimy;

		params = _params;

		cur = new TimeLayer2D(grid->dimx, grid->dimy, (FTYPE)grid->dx, (FTYPE)grid->dy);
		next = new TimeLayer2D(grid->dimx, grid->dimy, (FTYPE)grid->dx, (FTYPE)grid->dy);
		
		temp = new TimeLayer2D(grid->dimx, grid->dimy, (FTYPE)grid->dx, (FTYPE)grid->dy);
		next_w = new TimeLayer2D(grid->dimx, grid->dimy, (FTYPE)grid->dx, (FTYPE)grid->dy);
		q[0] = new ScalarField2D(grid->dimx, grid->dimy, (FTYPE)grid->dx, (FTYPE)grid->dy);
		q[1] = new ScalarField2D(grid->dimx, grid->dimy, (FTYPE)grid->dx, (FTYPE)grid->dy);
		div = new ScalarField2D(grid->dimx, grid->dimy, (FTYPE)grid->dx, (FTYPE)grid->dy);

		int maxInner = dimx * dimy;
		int maxBound = dimx * dimy;
		innerIndices = new int[maxInner * 2];
		boundIndices = new int[maxBound * 2];

		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
			{
				cur->U(i, j) = grid->GetData(i, j).vel.x;
				cur->V(i, j) = grid->GetData(i, j).vel.y;
				cur->T(i, j) = grid->GetData(i, j).T;
			}
		
		cur->CopyAllto(grid, next);
		cur->CopyAllto(grid, temp);
	}

	void StableSolver2D::SolveU(FTYPE dt, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next)
	{
		// eval new layer
		for (int ind = 0; ind < numInner; ind++)
		{
			int i = innerIndices[ind * 2 + 0];
			int j = innerIndices[ind * 2 + 1];
			next->U(i, j) = cur->U(i, j) + dt * ( 
				- temp->U(i, j) * temp->Ux(i, j) 
				- temp->V(i, j) * temp->Uy(i, j) 
				+ params.v_vis * (temp->Uxx(i, j) + temp->Uyy(i, j)) );
		}
	}

	void StableSolver2D::SolveV(FTYPE dt, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next)
	{
		// eval new layer
		for (int ind = 0; ind < numInner; ind++)
		{
			int i = innerIndices[ind * 2 + 0];
			int j = innerIndices[ind * 2 + 1];
			next->V(i, j) = cur->V(i, j) + dt * ( 
				- temp->U(i, j) * temp->Vx(i, j) 
				- temp->V(i, j) * temp->Vy(i, j) 
				+ params.v_vis * (temp->Vxx(i, j) + temp->Vyy(i, j)) );
		}
	}

	void StableSolver2D::Project(FTYPE dt, TimeLayer2D *w, TimeLayer2D *proj)
	{
		div->ClearZero();
		
		// eval divergence field first (optimization)
		#pragma omp parallel default(none) firstprivate(w)
		{
			#pragma omp for
			for (int ind = 0; ind < numInner; ind++)
			{
				int i = innerIndices[ind * 2 + 0];
				int j = innerIndices[ind * 2 + 1];
				div->U(i, j) = w->Ux(i, j) + w->Vy(i, j);
			}
		}

		// need to solve poisson equation : 
		//		laplas(q) = q_xx + q_yy = div(w)
		// simple discetization:
		//		(q[i+1,j] - 2q*[i,j] + q[i-1,j])/dx2 + (q[i,j+1] - 2q*[i,j] + q[i,j-1])/dy2 = div
		
		double dx2 = grid->dx * grid->dx;
		double dy2 = grid->dy * grid->dy;
		double rcp_dxdy2 = 0.5 / (dx2 + dy2);

		double err, i0, i1, j0, j1;
		int cur = 0;	
		q[cur]->ClearZero();
		do
		{
			err = 0.0;			
	
			// boundary cells 
			for (int ind = 0; ind < numBound; ind++)
			{
				int i = boundIndices[ind * 2 + 0];
				int j = boundIndices[ind * 2 + 1];
			
				// newman conditions
				if (grid->GetType(i-1, j) == NODE_IN) i0 = q[cur]->U(i-1, j); else i0 = q[cur]->U(i+1, j);
				if (grid->GetType(i+1, j) == NODE_IN) i1 = q[cur]->U(i+1, j); else i1 = q[cur]->U(i-1, j);
				if (grid->GetType(i, j-1) == NODE_IN) j0 = q[cur]->U(i, j-1); else j0 = q[cur]->U(i, j+1);
				if (grid->GetType(i, j+1) == NODE_IN) j1 = q[cur]->U(i, j+1); else j1 = q[cur]->U(i, j-1);
	
				double q_new = rcp_dxdy2 * ((i0 + i1) * dy2 + (j0 + j1) * dx2 - div->U(i, j) * dx2 * dy2);		
				double cur_err = abs((q_new - q[cur]->U(i, j)) / q_new);
				err = max(cur_err, err);
				q[cur]->U(i, j) = (FTYPE)q_new;
			}

			// inner cells
			for (int ind = 0; ind < numInner; ind++)
			{
				int i = innerIndices[ind * 2 + 0];
				int j = innerIndices[ind * 2 + 1];
			
				i0 = q[cur]->U(i-1, j); 
				i1 = q[cur]->U(i+1, j); 
				j0 = q[cur]->U(i, j-1); 
				j1 = q[cur]->U(i, j+1); 

				double q_new = rcp_dxdy2 * ((i0 + i1) * dy2 + (j0 + j1) * dx2 - div->U(i, j) * dx2 * dy2);		
				double cur_err = abs((q_new - q[cur]->U(i, j)) / q_new);
				err = max(cur_err, err);
				q[cur]->U(i, j) = (FTYPE)q_new;
			}

		} while (err >= POISSON_ERR_THRESHOLD);	// exit only when the error is small enough
		
		// force divergence free field: 
		//		proj = w - grad(q)

		for (int ind = 0; ind < numInner; ind++)
		{
			int i = innerIndices[ind * 2 + 0];
			int j = innerIndices[ind * 2 + 1];
			Vec2D g = q[cur]->grad(i, j);
			proj->U(i, j) = w->U(i, j) - g.x;
			proj->V(i, j) = w->V(i, j) - g.y;
		}
	}

	void StableSolver2D::PrepareIndices()
	{
		numInner = 0;
		numBound = 0;
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
			{
				switch (grid->GetType(i, j))
				{
				case NODE_IN:
					innerIndices[numInner * 2 + 0] = i;
					innerIndices[numInner * 2 + 1] = j;
					numInner++; 
					break;
				case NODE_BOUND:
				case NODE_VALVE:
					boundIndices[numBound * 2 + 0] = i;
					boundIndices[numBound * 2 + 1] = j;
					numBound++; 
					break;
				}
			}
	}

	void StableSolver2D::TimeStep(FTYPE dt, int num_global, int num_local)
	{
		cur->CopyAllto(grid, temp);

		PrepareIndices();

		// do global iterations
		int it;
		double err = next->EvalDivError(grid);

		for (it = 0; (it < num_global) || (err > DIV_ERR_THRESHOLD); it++)
		{
			cur->CopyAllto(grid, next_w);

			// solve equations
			SolveU(dt, cur, temp, next_w);
			SolveV(dt, cur, temp, next_w);
			Project(dt, next_w, next);

			err = next->EvalDivError(grid);
			
			// update non-linear parameters
			next->MergeAllto(grid, temp, NODE_IN);
			
			if (it > MAX_GLOBAL_ITERS) 
			{
				printf("\nExceeded max number of iterations (%i)\n", MAX_GLOBAL_ITERS); 
				exit(1); 
			}

			if (err > DIV_ERR_THRESHOLD * 10)
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

	void StableSolver2D::FreeMemory()
	{
		if (cur != NULL) delete cur;
		if (next != NULL) delete next;
		
		if (temp != NULL) delete temp;
		if (next_w != NULL) delete next_w;
		if (q[0] != NULL) delete q[0];
		if (q[1] != NULL) delete q[1];
		if (div != NULL) delete div;

		if (innerIndices != NULL) delete [] innerIndices;
		if (boundIndices != NULL) delete [] boundIndices;
	}

	StableSolver2D::StableSolver2D()
	{
		grid = NULL;

		cur = NULL;
		next = NULL;
		
		temp = NULL;
		next_w = NULL;
		q[0] = NULL;
		q[1] = NULL;
		div = NULL;

		innerIndices = NULL;
		boundIndices = NULL;
	}

	StableSolver2D::~StableSolver2D()
	{
		FreeMemory();
	}
}