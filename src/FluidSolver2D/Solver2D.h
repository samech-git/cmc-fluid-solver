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

#pragma once

#include "Grid2D.h"
#include "TimeLayer2D.h"

namespace FluidSolver2D
{
	class Solver2D
	{
	public:
		virtual void Init(Grid2D* grid, FluidParams &params) = 0;
		virtual void TimeStep(FTYPE dt, int num_global, int num_local) = 0;
		
		void GetLayer(Vec2D *v, double *T, int outdimx = 0, int outdimy = 0);
		void SetLayer(Vec2D *v, double *T);

		void UpdateBoundaries();
		void SetGridBoundaries();
		void ClearOutterCells();

		Grid2D *grid;

	protected:
		int dimx, dimy;
		FluidParams params;
		TimeLayer2D *cur, *next;

		double EvalDivError(TimeLayer2D *cur);
	};
}