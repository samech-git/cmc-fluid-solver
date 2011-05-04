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

#include "..\Common\Geometry.h"
#include "..\Common\Algorithms.h"

#include "Solver2D.h"

#define ERR_THRESHOLD		0.1
#define MAX_GLOBAL_ITERS	100

namespace FluidSolver2D
{
	enum VarType { type_U, type_V, type_T };
	enum DirType { X, Y };

	struct Segment2D
	{
		int posx, posy;
		int endx, endy;
		int size;
		DirType dir; 
	};

	class AdiSolver2D : public Solver2D
	{
	public:
		void Init(Grid2D* _grid, FluidParams &_params);
		void TimeStep(FTYPE dt, int num_global, int num_local);

		AdiSolver2D();
		~AdiSolver2D();
	
	private:
		vector<Segment2D> listX, listY;

		TimeLayer2D *temp, *half, *next_local;
		FTYPE *a, *b, *c, *d, *x;

		void BuildMatrix(FTYPE dt, int i, int j, VarType var, DirType dir, FTYPE *a, FTYPE *b, FTYPE *c, FTYPE *d, int n, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *temp_local);
		void ApplyBC0(int i, int j, VarType var, FTYPE &b0, FTYPE &c0, FTYPE &d0);
		void ApplyBC1(int i, int j, VarType var, FTYPE &a1, FTYPE &b1, FTYPE &d1);
		
		void CreateSegments();	
		void SolveSegment(FTYPE dt, Segment2D seg, VarType var, DirType dir, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *temp_local, TimeLayer2D *next_local);
		void UpdateSegment(FTYPE *x, Segment2D seg, VarType var, TimeLayer2D *layer);
		
		void SolveDirection(FTYPE dt, int num_local, vector<Segment2D> &list, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next);

		void FreeMemory();
	};
}