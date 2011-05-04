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

#include "Solver3D.h"

namespace FluidSolver3D
{
	void Solver3D::GetLayer(Vec3D *v, double *T, int outdimx, int outdimy, int outdimz)
	{
		next->Clear(grid, NODE_OUT, MISSING_VALUE, MISSING_VALUE, MISSING_VALUE, MISSING_VALUE);
		next->FilterToArrays(v, T, outdimx, outdimy, outdimz);
	}

	void Solver3D::UpdateBoundaries()
	{
		cur->CopyFromGrid(grid, NODE_BOUND);
		cur->CopyFromGrid(grid, NODE_VALVE);
		cur->CopyLayerTo(grid, next, NODE_BOUND);
		cur->CopyLayerTo(grid, next, NODE_VALVE);
	}

	void Solver3D::SetGridBoundaries()
	{
		cur->CopyToGrid(grid);
	}

	void Solver3D::ClearOutterCells()
	{
		next->Clear(grid, NODE_OUT, 0.0, 0.0, 0.0, (FTYPE)grid->baseT);
	}
}