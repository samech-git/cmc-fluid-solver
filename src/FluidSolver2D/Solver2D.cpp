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

#include "Solver2D.h"

namespace FluidSolver2D
{
	void Solver2D::GetLayer(Vec2D *v, double *T, int outdimx, int outdimy)
	{
		if (outdimx == 0) outdimx = dimx;
		if (outdimy == 0) outdimy = dimy;

		for (int i = 0; i < outdimx; i++)
			for (int j = 0; j < outdimy; j++)
			{
				int x = (i * dimx / outdimx);
				int y = (j * dimy / outdimy);
				v[i * outdimy + j].x = next->U(x, y); 
				v[i * outdimy + j].y = next->V(x, y);
				T[i * outdimy + j] = next->T(x, y);
			}
	}

	void Solver2D::SetLayer(Vec2D *v, double *T)
	{
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				{
					cur->U(i, j) = v[i * dimy + j].x;
					cur->V(i, j) = v[i * dimy + j].y;
					cur->T(i, j) = (FTYPE)T[i * dimy + j];
				}		
	}

	void Solver2D::UpdateBoundaries()
	{
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				switch(grid->GetType(i, j))
				{
				case NODE_BOUND:
				case NODE_VALVE:
					cur->U(i, j) = grid->GetData(i, j).vel.x;
					cur->V(i, j) = grid->GetData(i, j).vel.y;
					cur->T(i, j) = grid->GetData(i, j).T;
					break;
				}
		cur->CopyAllto(grid, next, NODE_BOUND);
		cur->CopyAllto(grid, next, NODE_VALVE);
	}

	void Solver2D::SetGridBoundaries()
	{
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
			{
				CondData2D d = grid->GetData(i, j);
				grid->SetFieldData(i, j, CondData2D(d.type, d.cell, Vec2D(cur->U(i, j), cur->V(i, j)), 0));
			}
	}

	void Solver2D::ClearOutterCells()
	{
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				if (grid->GetType(i, j) == NODE_OUT)
				{
					next->U(i, j) = 0.0;
					next->V(i, j) = 0.0;
					next->T(i, j) = (FTYPE)grid->startT;
				}
	}
}