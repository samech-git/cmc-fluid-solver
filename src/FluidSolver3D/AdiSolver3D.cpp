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

#include "AdiSolver3D.h"

#ifdef _WIN32
#include "..\Common\Algorithms.h"
#elif __unix__
#include "../Common/Algorithms.h"
#endif

//#include "../Common/test_util.h"

namespace FluidSolver3D
{

	double  AdiSolver3D::sum_layer(char ch)
	{
		PARAplan *pplan = PARAplan::Instance();

		TimeLayer3D* layer;
		switch (ch)
		{
		case 'c':
			layer = cur;
			break;
		case 't':
			layer = temp;
			break;
		case 'n':
			layer = next;
			break;
		case 'h':
			layer = half;
			break;
		}

		double sum = 0.;
		//int ndimx = pplan->getLength1D();
		//sum += TestUtil::sumEllementsMultiGPU(layer->U->getMultiArray(), ndimx*dimy*dimz, layer->haloSize);
		//sum += TestUtil::sumEllementsMultiGPU(layer->V->getMultiArray(), ndimx*dimy*dimz, layer->haloSize);
		//sum += TestUtil::sumEllementsMultiGPU(layer->W->getMultiArray(), ndimx*dimy*dimz, layer->haloSize);
		//sum += TestUtil::sumEllementsMultiGPU(layer->T->getMultiArray(), ndimx*dimy*dimz, layer->haloSize);
		return sum;
	}


	void AdiSolver3D::debug(bool _debug)
	{
		ifdebug = _debug;
	}

	AdiSolver3D::AdiSolver3D()
	{
		grid = NULL;

		ifdebug = false;

		cur = NULL;
		temp = NULL;
		half = NULL;
		next = NULL;

		curT = NULL;
		tempT = NULL;
		nextT = NULL;

		a = NULL;
		b = NULL;
		c = NULL;
		d = NULL;
		x = NULL;

		d_a = NULL;
		d_b = NULL;
		d_c = NULL;
		d_d = NULL;
		d_x = NULL;

		h_listX = NULL;
		h_listY = NULL;
		h_listZ = NULL;

		d_listX = NULL;
		d_listY = NULL;
		d_listZ = NULL;

		d_node_listX = NULL;
		d_node_listY = NULL;
		d_node_listZ = NULL;
	}

	void AdiSolver3D::FreeMemory()
	{
		if (cur != NULL) delete cur;
		if (temp != NULL) delete temp;
		if (half != NULL) delete half;
		if (next != NULL) delete next;

		if (curT != NULL) delete curT;
		if (tempT != NULL) delete tempT;
		if (nextT != NULL) delete nextT;

		if (a != NULL) delete [] a;
		if (b != NULL) delete [] b;
		if (c != NULL) delete [] c;
		if (d != NULL) delete [] d;
		if (x != NULL) delete [] x;

		if (h_listX != NULL) delete [] h_listX;
		if (h_listY != NULL) delete [] h_listY;
		if (h_listZ != NULL) delete [] h_listZ;

		if (d_a != NULL) multiDevFree<FTYPE>(d_a);
		if (d_b != NULL) multiDevFree<FTYPE>(d_b);
		if (d_c != NULL) multiDevFree<FTYPE>(d_c);
		if (d_d != NULL) multiDevFree<FTYPE>(d_d);
		if (d_x != NULL) multiDevFree<FTYPE>(d_x);

		if (d_listX != NULL) multiDevFree<Segment3D>(d_listX);
		if (d_listY != NULL) multiDevFree<Segment3D>(d_listY);
		if (d_listZ != NULL) multiDevFree<Segment3D>(d_listZ);

		if (d_node_listX != NULL) multiDevFree<NodesBoundary3D>(d_node_listX);
		if (d_node_listY != NULL) multiDevFree<NodesBoundary3D>(d_node_listY);
		if (d_node_listZ != NULL) multiDevFree<NodesBoundary3D>(d_node_listZ);

		if (mpi_buf != NULL) gpuSafeCall( cudaFreeHost(mpi_buf), "cudaFreeHost");

		for (int i = 0; i < 3; i++ )
			delete [] numSegsGPU[i];
	}

	AdiSolver3D::~AdiSolver3D()
	{
		prof.PrintTimings(csvFormat);

		FreeMemory();
	}

	void AdiSolver3D::SetOptionsGPU(bool _transposeOpt, bool _decomposeOpt)
	{
		transposeOpt = _transposeOpt;
		decomposeOpt = _decomposeOpt;
	}

	void AdiSolver3D::Init(BackendType _backend, bool _csv, Grid3D* _grid, FluidParams &_params)
	{
		grid = _grid;
		backend = _backend;
		csvFormat = _csv;

		PARAplan* pplan = PARAplan::Instance();
		int dimxNodeOffset = pplan->getOffset1D();
		int dimxNode = pplan->getLength1D();

		dimx = grid->dimx;
		dimy = grid->dimy;
		dimz = grid->dimz;

		params = _params;

		int n = max(dimx, max(dimy, dimz));

		mpi_buf = NULL;

		h_listX = new Segment3D[grid->dimy * grid->dimz * MAX_SEGS_PER_ROW];
		h_listY = new Segment3D[grid->dimx * grid->dimz * MAX_SEGS_PER_ROW];
		h_listZ = new Segment3D[grid->dimx * grid->dimy * MAX_SEGS_PER_ROW];

		if( backend == GPU )
		{
			GPUplan *pGPUplan = GPUplan::Instance();
			for (int i = 0; i < 3; i++ )
				numSegsGPU[i] = new int[pplan->gpuNum()];
			// create GPU matrices
			multiDevAlloc<FTYPE>( d_a, dimxNode * n * n * MAX_SEGS_PER_ROW, true, 2 * n * n * MAX_SEGS_PER_ROW); // with halos for intercommunication );
			multiDevAlloc<FTYPE>( d_b, dimxNode * n * n * MAX_SEGS_PER_ROW, true, 2 * n * n * MAX_SEGS_PER_ROW);
			multiDevAlloc<FTYPE>( d_c, dimxNode * n * n * MAX_SEGS_PER_ROW, true, 2 * n * n * MAX_SEGS_PER_ROW);
			multiDevAlloc<FTYPE>( d_d, dimxNode * n * n * MAX_SEGS_PER_ROW, true, 2 * n * n * MAX_SEGS_PER_ROW);
			multiDevAlloc<FTYPE>( d_x, dimxNode * n * n * MAX_SEGS_PER_ROW, true, 2 * n * n * MAX_SEGS_PER_ROW);

			multiDevAlloc<Segment3D>( d_listX, grid->dimy * grid->dimz * MAX_SEGS_PER_ROW, false, 0 ); // each device has relevant list of segments
			multiDevAlloc<Segment3D>( d_listY, grid->dimx * grid->dimz * MAX_SEGS_PER_ROW, false, 0 );
			multiDevAlloc<Segment3D>( d_listZ, grid->dimx * grid->dimy * MAX_SEGS_PER_ROW, false, 0 );	

			multiDevAlloc<NodesBoundary3D>( d_node_listX, grid->dimy * grid->dimz * MAX_SEGS_PER_ROW, false, 0 ); // each device has relevant list of node bounds
			multiDevAlloc<NodesBoundary3D>( d_node_listY, grid->dimx * grid->dimz * MAX_SEGS_PER_ROW, false, 0 );
			multiDevAlloc<NodesBoundary3D>( d_node_listZ, grid->dimx * grid->dimy * MAX_SEGS_PER_ROW, false, 0 );		

#ifdef __PARA
		gpuSafeCall(cudaHostAlloc(&mpi_buf, sizeof(FTYPE) * 2 * n * n * MAX_SEGS_PER_ROW, cudaHostAllocDefault), "cudaHostAlloc");
		//mpi_buf = new FTYPE[2 * n * n * MAX_SEGS_PER_ROW];	
#endif

		}
		else
		{
			a = new FTYPE[n * n * n * MAX_SEGS_PER_ROW];
			b = new FTYPE[n * n * n * MAX_SEGS_PER_ROW];
			c = new FTYPE[n * n * n * MAX_SEGS_PER_ROW];
			d = new FTYPE[n * n * n * MAX_SEGS_PER_ROW];
			x = new FTYPE[n * n * n * MAX_SEGS_PER_ROW];

			transposeOpt = false;
			decomposeOpt = false;
		}
		int haloSize = 0; // no MPI support in CPU case
		if (backend == GPU)
			haloSize = grid->dimy * grid->dimz;
		cur = new TimeLayer3D(backend, grid, grid->dimy * grid->dimz);  //create slices with halos
		if (!transposeOpt)
			half = new TimeLayer3D(backend, dimxNode, grid->dimy, grid->dimz, (FTYPE)grid->dx, (FTYPE)grid->dy, (FTYPE)grid->dz, haloSize);
		next = new TimeLayer3D(backend, dimxNode, grid->dimy, grid->dimz, (FTYPE)grid->dx, (FTYPE)grid->dy, (FTYPE)grid->dz, haloSize);
		temp = new TimeLayer3D(backend, dimxNode, grid->dimy, grid->dimz, (FTYPE)grid->dx, (FTYPE)grid->dy, (FTYPE)grid->dz, haloSize);

#if (TRANSPOSE_OPT == 1)
		if (transposeOpt)
		{
			curT = new TimeLayer3D(backend, dimxNode, grid->dimz, grid->dimy, (FTYPE)grid->dx, (FTYPE)grid->dz, (FTYPE)grid->dy, haloSize);
			tempT = new TimeLayer3D(backend, dimxNode, grid->dimz, grid->dimy, (FTYPE)grid->dx, (FTYPE)grid->dz, (FTYPE)grid->dy, haloSize);
			nextT = new TimeLayer3D(backend, dimxNode, grid->dimz, grid->dimy, (FTYPE)grid->dx, (FTYPE)grid->dz, (FTYPE)grid->dy, haloSize);
		}
#endif
	}

	void AdiSolver3D::OutputSegmentsInfo(int num, Segment3D *list, char *filename)
	{
		PARAplan *pplan = PARAplan::Instance();
		if (pplan->rank() == 0)
		{
			FILE *file = NULL;
			fopen_s(&file, filename, "w");
			int num_elem = 0;
			for (int i = 0; i < num; i++) num_elem += list[i].size;
			fprintf(file, "num_systems = %i\nunknown_elements = %i\nall_elements = %i\n", num, num_elem-num*2, num_elem);
			for (int i = 0; i < num; i++)
				fprintf(file, "%i ", list[i].size);
			fclose(file);
		}
	}

	void AdiSolver3D::UpdateBoundaries()
	{
		switch (cur->hw)
		{
		case CPU:
			CreateSegments(); //put here for consistency - could be placed just before the solver
			cur->CopyFromGrid(grid, NODE_BOUND);
			cur->CopyFromGrid(grid, NODE_VALVE);
			break;
		case GPU:
			CreateSegments(); //required
			cur->CopyGridBoundary(Z, numSegsGPU[Z], d_listZ, d_node_listZ);
			cur->CopyGridBoundary(Y, numSegsGPU[Y], d_listY, d_node_listY);
			cur->CopyGridBoundary(X, numSegsGPU[X], d_listX, d_node_listX);
			break;
		}

		cur->CopyLayerTo(grid, next, NODE_BOUND);
		cur->CopyLayerTo(grid, next, NODE_VALVE);
	}

	void AdiSolver3D::TimeStep(FTYPE dt, int num_global, int num_local)
	{	
		PARAplan *pplan = PARAplan::Instance();
		//CreateSegments();	

		//output segments info to file for benchmarking tridiagonal solvers
		//OutputSegmentsInfo(numSegs[X], h_listX, "segsX.txt");
		//OutputSegmentsInfo(numSegs[Y], h_listY, "segsY.txt");
		//OutputSegmentsInfo(numSegs[Z], h_listZ, "segsZ.txt");

		// setup non-linear layer
		prof.StartEvent();
		cur->CopyLayerTo(temp);
		prof.StopEvent("CopyLayer");

		// create transposed cur array if opt is enabled
		if (transposeOpt)
		{
			prof.StartEvent();
			cur->Transpose(curT);
			next->Transpose(nextT); // UVA
			prof.StopEvent("Transpose");
		}
		// do global iterations		
		for (int it = 0; it < num_global; it++)
		{
			// alternating directions
			SolveDirection(Z, dt, num_local, h_listZ, d_listZ, d_node_listZ, cur, temp, next);
			if (transposeOpt)
			{
				SolveDirection(Y, dt, num_local, h_listY, d_listY, d_node_listY, next, temp, cur);
				SolveDirection(X, dt, num_local, h_listX, d_listX, d_node_listX, cur, temp, next);
			}
			else
			{
				SolveDirection(Y, dt, num_local, h_listY, d_listY, d_node_listY, next, temp, half);
				SolveDirection(X, dt, num_local, h_listX, d_listX, d_node_listX, half, temp, next);
			}

			 //update non-linear layer
			switch (backend)
			{
			case CPU:
#if !INTERNAL_MERGE_ENABLE
			case GPU:
#endif
				prof.StartEvent();
				next->MergeLayerTo(grid, temp, NODE_IN);
				prof.StopEvent("MergeLayer");
				break;
			}
		}
		// smooth results
		//temp->Smooth(grid, next, NODE_IN);

		// compute error
		prof.StartEvent();
		double err = next->EvalDivError(grid);
		prof.StopEvent("EvalDivError");

		// check & output error
		if (err > ERR_THRESHOLD) {
			printf("\nError is too big! %f\n", err);
			throw runtime_error("");
		}
		else
			if (pplan->rank() == 0)
			{
				printf("\rerr = %.8f,", err);
				fflush(stdout);
			}

		// clear cells for dynamic grid update
		prof.StartEvent();
		ClearOutterCells();
		prof.StopEvent("ClearLayer");

		// swap current/next pointers 
		TimeLayer3D *tmp = next;
		next = cur;
		cur = tmp;
	}

	template<DirType dir>
	void AdiSolver3D::CreateListSegments(int &numSeg, Segment3D *h_list, Segment3D **d_list, NodesBoundary3D **d_node_list, int dim1, int dim2, int dim3)
	{	
		PARAplan *pplan = PARAplan::Instance();
		int dimxNode = pplan->getLength1D();
		int dimxNodeOffset = pplan->getOffset1D();

		Segment3D *hh_list;
		NodesBoundary3D *h_node_list;
		NodesBoundary3D *hh_node_list;

		switch (dir)
		{
		case X:
			hh_list = new Segment3D[grid->dimy * grid->dimz * MAX_SEGS_PER_ROW];
			h_node_list = new NodesBoundary3D[grid->dimy * grid->dimz * MAX_SEGS_PER_ROW];
			hh_node_list = new NodesBoundary3D[grid->dimy * grid->dimz * MAX_SEGS_PER_ROW];
			break;
		case Y:
			hh_list = new Segment3D[grid->dimx * grid->dimz * MAX_SEGS_PER_ROW];
			h_node_list = new NodesBoundary3D[grid->dimx * grid->dimz * MAX_SEGS_PER_ROW];
			hh_node_list = new NodesBoundary3D[grid->dimx * grid->dimz * MAX_SEGS_PER_ROW];
			break;
		case Z:
			hh_list = new Segment3D[grid->dimx * grid->dimy * MAX_SEGS_PER_ROW];
			h_node_list = new NodesBoundary3D[grid->dimx * grid->dimy * MAX_SEGS_PER_ROW];
			hh_node_list = new NodesBoundary3D[grid->dimx * grid->dimy * MAX_SEGS_PER_ROW];
			break;
		}

		grid->GenerateListSegments(numSeg, hh_list, dim1, dim2, dim3, dir);	

		if(backend == GPU)
			grid->GenerateGridBoundaries(hh_node_list, numSeg, hh_list, transposeOpt);

		numSeg = _nodeSplitListSegments<dir>(h_list, h_node_list, numSeg, hh_list, hh_node_list, dimxNode, dimxNodeOffset);
		
		if (ifdebug)
		{
			printf("CreateListSegments: Node(%d), total number of segments per direction %d: %d\n", pplan->rank(), dir, numSeg);
			fflush(stdout);
		}

		if( backend == GPU )
		{
			// copy segments info to GPU as well
			GPUplan *pGPUplan = GPUplan::Instance();

			dimxNodeOffset = 0;
			for (int i = 0; i < pGPUplan->size(); i++)
			{
				dimxNode = pGPUplan->node(i)->getLength1D();
				int nSegGPU = _nodeSplitListSegments<dir>(hh_list, hh_node_list, numSeg, h_list, h_node_list, dimxNode, dimxNodeOffset);

				numSegsGPU[dir][i] = nSegGPU;
				dimxNodeOffset += dimxNode;
				
				if (ifdebug)
				{
					printf("CreateListSegments: Node(%d), number of segments per direction %d, device %d:  %d\n",pplan->rank(), dir, i, nSegGPU);
					fflush(stdout);
				}

				pGPUplan->setDevice(i); // ???
				cudaMemcpy( d_list[i], hh_list, sizeof(Segment3D) * nSegGPU, cudaMemcpyHostToDevice );
				cudaMemcpy( d_node_list[i], hh_node_list, sizeof(NodesBoundary3D) * nSegGPU, cudaMemcpyHostToDevice );
			}
		}
		delete [] hh_node_list;
		delete [] h_node_list;
		delete [] hh_list;
	}

template<DirType dir>
	int AdiSolver3D::_nodeSplitListSegments(Segment3D *dest_list, NodesBoundary3D *dest_node_list, int numSeg, Segment3D *src_list, NodesBoundary3D *src_node_list, int length, int offset)
	{
		Segment3D new_seg;
		NodesBoundary3D new_node_bound;
		int nSeg = 0;
		for (int i = 0; i < numSeg; i++)
		{ 
			new_seg = src_list[i];
			new_node_bound = src_node_list[i];
			new_seg.skipX = true;
			if ( new_seg.posx < offset + length && // check if segment falls into given device
						new_seg.endx >= offset )   //for non X segments endx=posx
			{						
				new_seg.skipX = false;
				if (new_seg.posx < offset)
					if (new_seg.endx >= offset + length)
						new_seg.type = UNBOUND;
					else
						if (new_seg.type == BOUND || new_seg.type == BOUND_END)
							new_seg.type = BOUND_END;
						else
							new_seg.type = UNBOUND;
				else
					if (new_seg.endx >= offset + length)
						if (new_seg.type == BOUND || new_seg.type == BOUND_START)
							new_seg.type = BOUND_START;
						else
							new_seg.type = UNBOUND;
				new_seg.posx = max(0, new_seg.posx  - offset); // Shift positions along X according to GPU node number
				new_seg.endx = min(length - 1, new_seg.endx - offset);
				new_seg.size = (new_seg.endx - new_seg.posx) + (new_seg.endy - new_seg.posy) + (new_seg.endz - new_seg.posz) + 1;

				dest_list[nSeg] = new_seg;
				dest_node_list[nSeg] = new_node_bound;
				nSeg++;
			}
			else
				switch (dir)
				{
					case X:
						new_seg.skipX = true; // we want to keep all segments along X for each node for segment id correspondence between different GPUs
						dest_list[nSeg] = new_seg;
						dest_node_list[nSeg] = new_node_bound;
						nSeg++;
						break;
				}
		}
		return nSeg;
	}

	void AdiSolver3D::CreateSegments()
	{
		prof.StartEvent();

		CreateListSegments<X>(numSegs[X], h_listX, d_listX, d_node_listX, dimx, dimy, dimz);
		CreateListSegments<Y>(numSegs[Y], h_listY, d_listY, d_node_listY, dimy, dimx, dimz);
		CreateListSegments<Z>(numSegs[Z], h_listZ, d_listZ, d_node_listZ, dimz, dimx, dimy);

		prof.StopEvent("CreateSegments");
	}

	void AdiSolver3D::SolveDirection(DirType dir, FTYPE dt, int num_local, Segment3D *h_list, Segment3D **d_list, NodesBoundary3D **d_node_list, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *next)
	{
		DirType dir_new = dir;
		TimeLayer3D *cur_new = cur;
		TimeLayer3D *temp_new = temp;
		TimeLayer3D *next_new = next;

		NodeType **d_node_types = grid->GetTypesGPU();

		// transpose non-linear layer if opt is enabled
		if( transposeOpt && dir == Z )
		{
			prof.StartEvent();			
			temp->Transpose(tempT);
			prof.StopEvent("Transpose");

			// set transposed direction
			dir_new = Z_as_Y;
			cur_new = curT;
			temp_new = tempT;
			next_new = nextT;
		}

		for (int it = 0; it < num_local; it++)
		{
			switch( backend )
			{
			case CPU:
				prof.StartEvent();
				#pragma omp parallel default(none) firstprivate(dt, dir) shared(h_list, cur, temp, next)
				{
					#pragma omp for
					for (int s = 0; s < numSegs[dir]; s++)
					{		
						SolveSegment(dt, s, h_list[s], type_U, dir, cur, temp, next);
						SolveSegment(dt, s, h_list[s], type_V, dir, cur, temp, next);
						SolveSegment(dt, s, h_list[s], type_W, dir, cur, temp, next);
						SolveSegment(dt, s, h_list[s], type_T, dir, cur, temp, next);			
					}
				}
				break;		
			case GPU:
				PARAplan *pplan = PARAplan::Instance();
				prof.StartEvent();
				temp_new->syncHalos(mpi_buf);
				switch( dir )
				{
				printf("Syncing Halos...\n");
				fflush(stdout);
				case X: prof.StopEvent("syncHalos_X"); break;
				case Y: prof.StopEvent("syncHalos_Y"); break;
				case Z: prof.StopEvent("syncHalos_Z"); break;
				}

				prof.StartEvent();

				SolveSegments_GPU(dt, params, numSegsGPU[dir], d_list, type_U, dir_new, d_node_list, d_node_types, cur_new, temp_new, next_new, d_a, d_b, d_c, d_d, d_x, decomposeOpt, numSegs[dir], mpi_buf);	
				SolveSegments_GPU(dt, params, numSegsGPU[dir], d_list, type_V, dir_new, d_node_list, d_node_types, cur_new, temp_new, next_new, d_a, d_b, d_c, d_d, d_x, decomposeOpt, numSegs[dir], mpi_buf);
				SolveSegments_GPU(dt, params, numSegsGPU[dir], d_list, type_W, dir_new, d_node_list, d_node_types, cur_new, temp_new, next_new, d_a, d_b, d_c, d_d, d_x, decomposeOpt, numSegs[dir], mpi_buf);
				SolveSegments_GPU(dt, params, numSegsGPU[dir], d_list, type_T, dir_new, d_node_list, d_node_types, cur_new, temp_new, next_new, d_a, d_b, d_c, d_d, d_x, decomposeOpt, numSegs[dir], mpi_buf);				
				break;
			}

			switch( dir )
			{
			case X: prof.StopEvent("SolveSegments_X"); break;
			case Y: prof.StopEvent("SolveSegments_Y"); break;
			case Z: prof.StopEvent("SolveSegments_Z"); break;
			}

			switch (backend)
			{
			case CPU:
#if !INTERNAL_MERGE_ENABLE
			case GPU:
#endif
				if( dir_new == Z_as_Y )
				{
					// update non-linear layer
					prof.StartEvent();
					nextT->MergeLayerTo(grid, tempT, NODE_IN, true);
					prof.StopEvent("MergeLayer");
				}
				else
				{
					// update non-linear layer
					prof.StartEvent();
					next->MergeLayerTo(grid, temp, NODE_IN);
					prof.StopEvent("MergeLayer");
				}
				break;
			}
		}

		// transpose temp and next layers to normal order
		if( dir_new == Z_as_Y )
		{
			prof.StartEvent();			
			nextT->Transpose(next);
			tempT->Transpose(temp);
			prof.StopEvent("Transpose");
		}
	}

void AdiSolver3D::SolveDirection_XY(FTYPE dt, int num_local, Segment3D *h_list, Segment3D **d_list, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *next)
{
}

	void AdiSolver3D::SolveSegment(FTYPE dt, int id, Segment3D seg, VarType var, DirType dir, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *next)
	{
		int n = seg.size;

		int max_n = max(dimx, max(dimy, dimz));
		FTYPE *a = this->a + id * max_n;
		FTYPE *b = this->b + id * max_n;
		FTYPE *c = this->c + id * max_n;
		FTYPE *d = this->d + id * max_n;
		FTYPE *x = this->x + id * max_n;

		ApplyBC0(seg.posx, seg.posy, seg.posz, var, b[0], c[0], d[0]);
		ApplyBC1(seg.endx, seg.endy, seg.endz, var, a[n-1], b[n-1], d[n-1]);
		BuildMatrix(dt, seg.posx, seg.posy, seg.posz, var, dir, a, b, c, d, n, cur, temp);
		
		SolveTridiagonal(a, b, c, d, x, n);
		
		UpdateSegment(x, seg, var, next);
	}

	void AdiSolver3D::UpdateSegment(FTYPE *x, Segment3D seg, VarType var, TimeLayer3D *layer)
	{
		int i = seg.posx;
		int j = seg.posy;
		int k = seg.posz;

		for (int t = 0; t < seg.size; t++)
		{
			switch (var)
			{
			case type_U: layer->U->elem(i, j, k) = x[t]; break;
			case type_V: layer->V->elem(i, j, k) = x[t]; break;
			case type_W: layer->W->elem(i, j, k) = x[t]; break;
			case type_T: layer->T->elem(i, j, k) = x[t]; break;
			}

			switch (seg.dir)
			{
			case X: i++; break;
			case Y: j++; break;
			case Z: k++; break;
			}
		}
	}

	void AdiSolver3D::BuildMatrix(FTYPE dt, int i, int j, int k, VarType var, DirType dir, FTYPE *a, FTYPE *b, FTYPE *c, FTYPE *d, int n, TimeLayer3D *cur, TimeLayer3D *temp)
	{
		FTYPE vis_dx2, vis_dy2, vis_dz2;
		FTYPE dx = cur->dx;
		FTYPE dy = cur->dy;
		FTYPE dz = cur->dz;

		switch (var)
		{
		case type_U:
		case type_V:
		case type_W:
			vis_dx2 = params.v_vis / (dx * dx);
			vis_dy2 = params.v_vis / (dy * dy);
			vis_dz2 = params.v_vis / (dz * dz);
			break;
		case type_T:
			vis_dx2 = params.t_vis / (dx * dx);
			vis_dy2 = params.t_vis / (dy * dy);
			vis_dz2 = params.t_vis / (dz * dz);
			break;
		}
		
		for (int p = 1; p < n-1; p++)
		{
			switch (dir)
			{
			case X:		
				a[p] = - temp->U->elem(i+p, j, k) / (2 * dx) - vis_dx2; 
				b[p] = 3 / dt  +  2 * vis_dx2; 
				c[p] = temp->U->elem(i+p, j, k) / (2 * dx) - vis_dx2; 
				
				switch (var)	
				{
				case type_U: d[p] = cur->U->elem(i+p, j, k) * 3 / dt - params.v_T * temp->T->d_x(i+p, j, k); break;
				case type_V: d[p] = cur->V->elem(i+p, j, k) * 3 / dt; break;
				case type_W: d[p] = cur->W->elem(i+p, j, k) * 3 / dt; break;
				case type_T: d[p] = cur->T->elem(i+p, j, k) * 3 / dt + params.t_phi * temp->DissFuncX(i+p, j, k); break;
				}	
				break;

			case Y:
				a[p] = - temp->V->elem(i, j+p, k) / (2 * dy) - vis_dy2; 
				b[p] = 3 / dt  +  2 * vis_dy2; 
				c[p] = temp->V->elem(i, j+p, k) / (2 * dy) - vis_dy2; 
				
				switch (var)	
				{
				case type_U: d[p] = cur->U->elem(i, j+p, k) * 3 / dt; break;
				case type_V: d[p] = cur->V->elem(i, j+p, k) * 3 / dt - params.v_T * temp->T->d_y(i, j+p, k); break;
				case type_W: d[p] = cur->W->elem(i, j+p, k) * 3 / dt; break;
				case type_T: d[p] = cur->T->elem(i, j+p, k) * 3 / dt + params.t_phi * temp->DissFuncY(i, j+p, k); break;
				}
				break;

			case Z:
				a[p] = - temp->W->elem(i, j, k+p) / (2 * dz) - vis_dz2; 
				b[p] = 3 / dt  +  2 * vis_dz2; 
				c[p] = temp->W->elem(i, j, k+p) / (2 * dz) - vis_dz2; 
				
				switch (var)	
				{
				case type_U: d[p] = cur->U->elem(i, j, k+p) * 3 / dt; break;
				case type_V: d[p] = cur->V->elem(i, j, k+p) * 3 / dt; break;
				case type_W: d[p] = cur->W->elem(i, j, k+p) * 3 / dt - params.v_T * temp->T->d_z(i, j, k+p); break;
				case type_T: d[p] = cur->T->elem(i, j, k+p) * 3 / dt + params.t_phi * temp->DissFuncZ(i, j, k+p); break;
				}
				break;
			}
		}
	}

	void AdiSolver3D::ApplyBC0(int i, int j, int k, VarType var, FTYPE &b0, FTYPE &c0, FTYPE &d0)
	{
		if ((var == type_T && grid->GetBC_temp(i, j, k) == BC_FREE) ||
			(var != type_T && grid->GetBC_vel(i, j, k) == BC_FREE))
		{
			// free: f(0) = 2 * f(1) - f(2)
			b0 = 2.0; 
			c0 = -1.0; 
			d0 = 0.0; 
		}
		else
		{
			// no-slip: f(0) = f(1)
			b0 = 1.0; 
			c0 = 0.0; 
			switch (var)
			{
			case type_U: d0 = (FTYPE)grid->GetVel(i, j, k).x; break;
			case type_V: d0 = (FTYPE)grid->GetVel(i, j, k).y; break;
			case type_W: d0 = (FTYPE)grid->GetVel(i, j, k).z; break;
			case type_T: d0 = (FTYPE)grid->GetT(i, j, k); break;
			}
		}
	}

	void AdiSolver3D::ApplyBC1(int i, int j, int k, VarType var, FTYPE &a1, FTYPE &b1, FTYPE &d1)
	{
		if ((var == type_T && grid->GetBC_temp(i, j, k) == BC_FREE) ||
			(var != type_T && grid->GetBC_vel(i, j, k) == BC_FREE))
		{
			// free: f(N) = 2 * f(N-1) - f(N-2)
			a1 = -1.0; 
			b1 = 2.0; 
			d1 = 0.0;
		}
		else
		{
			// no-slip: f(N) = f(N-1)
			a1 = 0.0; 
			b1 = 1.0; 
			switch (var)
			{
			case type_U: d1 = (FTYPE)grid->GetVel(i, j, k).x; break;
			case type_V: d1 = (FTYPE)grid->GetVel(i, j, k).y; break;
			case type_W: d1 = (FTYPE)grid->GetVel(i, j, k).z; break;
			case type_T: d1 = (FTYPE)grid->GetT(i, j, k); break;
			}
		}
	}	
}
