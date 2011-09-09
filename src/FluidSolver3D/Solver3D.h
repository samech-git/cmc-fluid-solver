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

#include "Grid3D.h"
#include "TimeLayer3D.h"

namespace FluidSolver3D
{
	class Solver3D
	{
	public:
		virtual void Init(BackendType backend, bool csvFormat, Grid3D* grid, FluidParams &params) = 0;
		virtual void TimeStep(FTYPE dt, int num_global, int num_local) = 0;
		virtual double sum_layer(char ch) = 0;
		virtual void debug(bool ifdebug) = 0;
		
		void GetLayer(Vec3D *v, double *T, int outdimx = 0, int outdimy = 0, int outdimz = 0);

		void UpdateBoundaries();
		void SetGridBoundaries();
		void ClearOutterCells();

		Grid3D *grid;

		virtual ~Solver3D() {};

	protected:
		int dimx, dimy, dimz;
		FluidParams params;
		TimeLayer3D *cur, *next;		

		double EvalDivError(TimeLayer3D *cur);
	};
}
