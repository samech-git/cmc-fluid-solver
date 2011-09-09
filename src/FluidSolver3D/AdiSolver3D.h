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

#define PROFILE_ENABLE		1
#define BLOCKING_SOLVER_ENABLE 1
#define INTERNAL_MERGE_ENABLE 0

#ifdef _WIN32
#include "..\Common\Profiler.h"
#elif __unix__
#include "../Common/Profiler.h"
#endif

#include "Solver3D.h"

#define ERR_THRESHOLD		0.01

using namespace Common;

namespace FluidSolver3D
{
	enum VarType { type_U, type_V, type_W, type_T };

	extern void SolveSegments_GPU( FTYPE dt, FluidParams params, int* num_seg, Segment3D **segs, VarType var, DirType dir, Node **nodes, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *next,
								   FTYPE **d_a, FTYPE **d_b, FTYPE **d_c, FTYPE **d_d, FTYPE **d_x, bool decomposeOpt, int numSegs, FTYPE *mpi_buf = NULL);

	class AdiSolver3D : public Solver3D
	{
	public:
		AdiSolver3D();
		~AdiSolver3D();

		void Init(BackendType backend, bool _csv, Grid3D* _grid, FluidParams &_params);
		void TimeStep(FTYPE dt, int num_global, int num_local);
		void SetOptionsGPU(bool _transposeOpt, bool _decomposeOpt);
		double sum_layer(char ch);
		void debug(bool ifdebug);

	private:
		Profiler prof;
		bool csvFormat;
		BackendType backend;

		// options, optimizations
		bool transposeOpt, decomposeOpt;
		
		int numSegs[3];
		int* numSegsGPU[3];  // segments per direction per GPU
		Segment3D *h_listX, *h_listY, *h_listZ;				// segments in CPU mem
		Segment3D **d_listX, **d_listY, **d_listZ;		// segments in multiple GPU mem 

		TimeLayer3D *temp, *half;
		TimeLayer3D *curT, *tempT, *nextT;			// for transpose GPU optimization

		FTYPE *mpi_buf;

		FTYPE *a, *b, *c, *d, *x;									// matrices in CPU mem
		FTYPE **d_a, **d_b, **d_c, **d_d, **d_x; // same matrices in GPU mem

		bool ifdebug; // flag for printing stats only;

		void BuildMatrix(FTYPE dt, int i, int j, int k, VarType var, DirType dir, FTYPE *a, FTYPE *b, FTYPE *c, FTYPE *d, int n, TimeLayer3D *cur, TimeLayer3D *temp);
		void ApplyBC0(int i, int j, int k, VarType var, FTYPE &b0, FTYPE &c0, FTYPE &d0);
		void ApplyBC1(int i, int j, int k, VarType var, FTYPE &a1, FTYPE &b1, FTYPE &d1);
		
		template<DirType dir>
		void CreateListSegments(int &numSeg, Segment3D *h_list, Segment3D **d_list, int dim1, int dim2, int dim3);
		template<DirType dir>
		int _nodeSplitListSegments(Segment3D *dest_list, int numSeg, Segment3D *src_list, int dimx, int dimxOffset);
		
		void OutputSegmentsInfo(int num, Segment3D *list, char *filename);
		void CreateSegments();

		void SolveSegment(FTYPE dt, int id, Segment3D seg, VarType var, DirType dir, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *next);
		void UpdateSegment(FTYPE *x, Segment3D seg, VarType var, TimeLayer3D *layer);
		
		void SolveDirection(DirType dir, FTYPE dt, int num_local, Segment3D *h_list, Segment3D **d_list, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *next);
		void SolveDirection_XY(FTYPE dt, int num_local, Segment3D *h_list, Segment3D **d_list, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *next);

		void FreeMemory();
	};
}
