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

#define SEG_BLOCK_DIM_X		32
#define SEG_BLOCK_DIM_Y		8

#if( __CUDA_ARCH__ < 120 )
#define SOLVER_BLOCK_DIM	128
#else
#define SOLVER_BLOCK_DIM	256
#endif

namespace FluidSolver3D
{
	struct FluidParamsGPU
	{
		FTYPE vis_dx2;
		FTYPE dt, dx, dy, dz;
		FTYPE v_T, t_phi;

		//FluidParamsGPU() :  dt(0.0), dx(0.0), dy(0.0), dz(0.0){;};
		FluidParamsGPU( VarType var, DirType dir, FTYPE _dt, FTYPE _dx, FTYPE _dy, FTYPE _dz, FluidParams _params ) :
			dt(_dt), dx(_dx), dy(_dy), dz(_dz)
		{
			switch (var)
			{
			case type_U:
			case type_V:
			case type_W:
				switch (dir)
				{
				case X:	vis_dx2 = _params.v_vis / (dx * dx); break;
				case Y: vis_dx2 = _params.v_vis / (dy * dy); break;
				case Z: case Z_as_Y: vis_dx2 = _params.v_vis / (dz * dz); break;
				}
				break;
			case type_T:
				switch (dir)
				{
				case X:	vis_dx2 = _params.t_vis / (dx * dx); break;
				case Y: vis_dx2 = _params.t_vis / (dy * dy); break;
				case Z: case Z_as_Y:  vis_dx2 = _params.t_vis / (dz * dz); break;
				}
				break;
			}
						
			v_T = _params.v_T;
			t_phi = _params.t_phi;
		}
	};

#if 1
	// interleave matrix arrays for better memory access
	#define get(a, elem_id)			a[id + (elem_id + 1) * max_n_max_n * MAX_SEGS_PER_ROW]
#else
	// sequential layout - bad access pattern,  currently not implemented for MGPU
//	#define get(a, elem_id)			a[elem_id + id * max_n]
#endif

template<int dir, int swipe, int var>
	__device__ void solve_tridiagonal(FluidParamsGPU params, FTYPE *c, FTYPE *d, FTYPE *x, Segment3D &seg, int num, NodesBoundary3D &nodesBound, TimeLayer3D_GPU &cur, TimeLayer3D_GPU &temp, int id, int num_seg, int max_n_max_n, int dimX, SegmentType segType = BOUND)
	{
		FTYPE a_val, b_val, c_val, d_val;
		int i = seg.posx;
		int j = seg.posy;
		int k = seg.posz;

		switch (swipe)
		{
			case ALL:
			case FORWARD:
			{
				// apply boundary 0 conditions:
				switch (dir)
				{
				case X:
					if (segType == UNBOUND || segType == BOUND_END)
						break;
				default:
					if ((var == type_T && nodesBound.first.bc_temp == BC_FREE) ||
						(var != type_T && nodesBound.first.bc_vel == BC_FREE))
					{
						// free: f(0) = 2 * f(1) - f(2)
						b_val = 2.0; 
						c_val = -1.0; 
						d_val = 0.0;
					}
					else
					{
						// no-slip: f(0) = f(1)
						b_val = 1.0; 
						c_val = 0.0; 

						switch (var)
						{
							case type_U: d_val = nodesBound.first.v.x; break;  
							case type_V: d_val = nodesBound.first.v.y; break;
							case type_W: d_val = nodesBound.first.v.z; break;
							case type_T: d_val = nodesBound.first.T; break;
						}
					}
				}

				int start = 1;
				int end = num - 1;

				switch (dir)
				{
				case X:
						switch (segType)
						{
						case UNBOUND:
							start = 0;
							end = num;
							break;
						case BOUND_END:
							start = 0;
							//get(c, num-1) = 0.0;
							break;
						case BOUND:
							//get(c, num-1) = 0.0;
						case BOUND_START:
							get(c,0) = c_val / b_val;
							get(d,0) = d_val / b_val;
							end = num;
							break;
						}
						break;
				default:
					get(c,0) = c_val / b_val;
					get(d,0) = d_val / b_val;
					//get(c,num-1) = 0.0;
					break;
				}

				for (int p = start; p < end; p++)
				{
					// Build matrix:
					switch (dir)
					{
					case X:		
						a_val = - temp.elem(temp.u, i+p, j, k) / (2 * params.dx) - params.vis_dx2; 
						b_val = 3 / params.dt  +  2 * params.vis_dx2; 
						c_val = temp.elem(temp.u, i+p, j, k) / (2 * params.dx) - params.vis_dx2; 
				
						switch (var)	
						{
						case type_U: d_val = cur.elem(cur.u, i+p, j, k) * 3 / params.dt - params.v_T * temp.d_x(temp.T, params.dx, i+p, j, k); break;
						case type_V: d_val = cur.elem(cur.v, i+p, j, k) * 3 / params.dt; break;
						case type_W: d_val = cur.elem(cur.w, i+p, j, k) * 3 / params.dt; break;
						case type_T: d_val = cur.elem(cur.T, i+p, j, k) * 3 / params.dt + params.t_phi * temp.DissFuncX(params.dx, params.dy, params.dz, i+p, j, k); break;
						}	
						break;

					case Y:
						a_val = - temp.elem(temp.v, i, j+p, k) / (2 * params.dy) - params.vis_dx2; 
						b_val = 3 / params.dt  +  2 * params.vis_dx2; 
						c_val = temp.elem(temp.v, i, j+p, k) / (2 * params.dy) - params.vis_dx2; 
				
						switch (var)	
						{
						case type_U: d_val = cur.elem(cur.u, i, j+p, k) * 3 / params.dt; break;
						case type_V: d_val = cur.elem(cur.v, i, j+p, k) * 3 / params.dt - params.v_T * temp.d_y(temp.T, params.dy, i, j+p, k); break;
						case type_W: d_val = cur.elem(cur.w, i, j+p, k) * 3 / params.dt; break;
						case type_T: d_val = cur.elem(cur.T, i, j+p, k) * 3 / params.dt + params.t_phi * temp.DissFuncY(params.dx, params.dy, params.dz, i, j+p, k); break;
						}
						break;

					case Z:
						a_val = - temp.elem(temp.w, i, j, k+p) / (2 * params.dz) - params.vis_dx2; 
						b_val = 3 / params.dt  +  2 * params.vis_dx2; 
						c_val = temp.elem(temp.w, i, j, k+p) / (2 * params.dz) - params.vis_dx2; 
				
						switch (var)	
						{
						case type_U: d_val = cur.elem(cur.u, i, j, k+p) * 3 / params.dt; break;
						case type_V: d_val = cur.elem(cur.v, i, j, k+p) * 3 / params.dt; break;
						case type_W: d_val = cur.elem(cur.w, i, j, k+p) * 3 / params.dt - params.v_T * temp.d_z(temp.T, params.dz, i, j, k+p); break;
						case type_T: d_val = cur.elem(cur.T, i, j, k+p) * 3 / params.dt + params.t_phi * temp.DissFuncZ(params.dx, params.dy, params.dz, i, j, k+p); break;
						}
						break;

					case Z_as_Y:
						a_val = - temp.elem(temp.w, i, k+p, j) / (2 * params.dz) - params.vis_dx2; 
						b_val = 3 / params.dt  +  2 * params.vis_dx2; 
						c_val = temp.elem(temp.w, i, k+p, j) / (2 * params.dz) - params.vis_dx2; 
				
						switch (var)	
						{
						case type_U: d_val = cur.elem(cur.u, i, k+p, j) * 3 / params.dt; break;
						case type_V: d_val = cur.elem(cur.v, i, k+p, j) * 3 / params.dt; break;
						case type_W: d_val = cur.elem(cur.w, i, k+p, j) * 3 / params.dt - params.v_T * temp.d_y(temp.T, params.dz, i, k+p, j); break;
						case type_T: d_val = cur.elem(cur.T, i, k+p, j) * 3 / params.dt + params.t_phi * temp.DissFuncZ_as_Y(params.dx, params.dz, params.dy, i, k+p, j); break;
						}
						break;
					} // end of build matrix

					// forward solver step:
					get(c,p) = c_val / (b_val - a_val * get(c,p-1));
					get(d,p) = (d_val - get(d,p-1) * a_val) / (b_val - a_val * get(c,p-1));

					//c_cur = c_val / (b_val - a_val * c_m1);
					//d_cur = (d_val - d_m1 * a_val) / (b_val - a_val * c_m1);
				}

				// apply boundary 1 conditions:
				switch (dir)
				{
				case X:
					if (segType == UNBOUND || segType == BOUND_START)
						break;
				default:
					if ((var == type_T && nodesBound.last.bc_temp == BC_FREE) ||
						(var != type_T && nodesBound.last.bc_vel == BC_FREE))
					{
						// free: f(N) = 2 * f(N-1) - f(N-2)
						a_val = -1.0;
						b_val = 2.0;
						d_val = 0.0;
					}
					else
					{
						// no-slip: f(N) = f(N-1)
						a_val = 0.0; 
						b_val = 1.0; 

						switch (var)
						{
							case type_U: d_val = nodesBound.last.v.x; break;
							case type_V: d_val = nodesBound.last.v.y; break;
							case type_W: d_val = nodesBound.last.v.z; break;
							case type_T: d_val = nodesBound.last.T; break;
						}
					}
					// move boundary setting in here
					get(c, num-1) = 0.0;
					get(d,num-1) = (d_val - get(d,num-2) * a_val) / (b_val - a_val * get(c,num-2));
				}

				switch (dir)
				{
				case X:
					// copy segments' ends to halo area
					get(c, dimX-1) = get(c, num-1);
					get(d, dimX-1) = get(d, num-1);
					break;
				}

				break; // SWIPE
			}			
		}

				switch (dir)
				{
				case Y:
					;
					break;
				}

		switch (swipe)
		{
			case ALL:
			case BACK:
			{
				int end = num - 1;

				switch (dir)
				{
				case X:
					switch (segType)
					{
					case UNBOUND:
						end = num;
						get(x,num) = get(x,dimX);
						break;
					case BOUND_START:
						end = num;
						get(x,num) = get(x,dimX);
						break;
					case BOUND_END:
						get(x, num-1) = get(d,num-1);
						break;
					case BOUND:
						get(x, num-1) = get(d,num-1);
						break;
					}
					break;
				default:
					get(x, num-1) = get(d,num-1);
					break;
				}

				for (int i = end-1; i >= 0; i--) 
					get(x,i) = get(d,i) -  get(c,i) * get(x, i+1);
				break;
			}
		}
	}

	template<int dir, int var>
	__device__ void update_segment( FTYPE *x, Segment3D &seg, NodeType *nodes, TimeLayer3D_GPU &temp, TimeLayer3D_GPU &layer, int id, int num_seg, int max_n_max_n )
	{
		int i = seg.posx;
		int j = seg.posy;
		int k = seg.posz;

		if( dir == Z_as_Y ) 
		{
			k = seg.posy;
			j = seg.posz;
		}

		for (int t = 0; t < seg.size; t++)
		{
			switch (var)
			{
			case type_U: layer.elem(layer.u, i, j, k) = get(x,t); break;
			case type_V: layer.elem(layer.v, i, j, k) = get(x,t); break;
			case type_W: layer.elem(layer.w, i, j, k) = get(x,t); break;
			case type_T: layer.elem(layer.T, i, j, k) = get(x,t); break;
			}
#if INTERNAL_MERGE_ENABLE == 1
			//if (t == 0)
			//	if (nodes.first.type != NODE_IN)
			//		continue;
			//if (t == seg.size-1)
			//	if (nodes.last.type != NODE_IN)
			//		continue;

			int idn;
			idn = i * layer.dimy * layer.dimz;
			switch(dir)
			{
			case Z_as_Y:
				idn += k * layer.dimy + j;
				break;
			default:
				idn += j * layer.dimz + k;
			}
			//int idn = i * layer.dimy * layer.dimz + j * layer.dimz + k;
			if (nodes[idn] == NODE_IN)
			switch (var)
			{
			case type_U: temp.elem(temp.u, i, j, k) = (temp.elem(temp.u, i, j, k) +  get(x,t) ) / 2; break;
			case type_V: temp.elem(temp.v, i, j, k) = (temp.elem(temp.v, i, j, k) +  get(x,t) ) / 2; break;
			case type_W: temp.elem(temp.w, i, j, k) = (temp.elem(temp.w, i, j, k) +  get(x,t) ) / 2; break;
			case type_T: temp.elem(temp.T, i, j, k) = (temp.elem(temp.T, i, j, k) +  get(x,t) ) / 2; break;
			}
#endif

			switch (dir)
			{
			case X: i++; break;
			case Y: j++; break;
			case Z: k++; break;
			case Z_as_Y: j++; break;
			}
		}
	}

	template<int dir, int swipe, int var>
	__global__ void solve_segments_var(FluidParamsGPU p, int num_seg, Segment3D *segs, NodesBoundary3D *nodesBounds, NodeType *nodeTypes, TimeLayer3D_GPU cur, TimeLayer3D_GPU temp, TimeLayer3D_GPU next,
									FTYPE *c, FTYPE *x, int  max_n_max_n, int dimX, int id_shift = 0 )
	{
		// fetch current segment info
		int id = id_shift + blockIdx.x * blockDim.x + threadIdx.x;
		if( id >= num_seg + id_shift) return;
		Segment3D &seg = segs[id];
		NodesBoundary3D &nodes = nodesBounds[id];

		int n = seg.size;

		switch (dir)
		{
		case X:
			if (seg.skipX)
				return;
		}

		solve_tridiagonal<dir, swipe, var>(p, c, x, x, seg, n, nodes, cur, temp, id, num_seg, max_n_max_n, dimX, seg.type);

		switch (swipe)
		{
		case ALL:
			update_segment<dir, var>(x, seg, nodeTypes, temp, next, id, num_seg, max_n_max_n);
		case BACK:			
			break;
		}
	}

	//template<DirType dir, int swipe>
	//void solve_segments(dim3 grid, dim3 block, FluidParamsGPU p, int num_seg, Segment3D *segs, NodesBoundary3D *nodesBounds, NodeType *nodeTypes, TimeLayer3D_GPU cur, TimeLayer3D_GPU temp, TimeLayer3D_GPU next,
	//								FTYPE *c, FTYPE *x, int  max_n_max_n, int dimX, int id_shift = 0)
	//{
	//	solve_segments<dir, swipe, type_U><<<grid, block>>>( *p[type_U], num_seg, segs, nodesBounds, nodeTypes, cur, temp, next, c, d_x[i] + typeID*varOffset, max_n_max_n, dimX );
	//}

	template<int dir, int var>
	__global__ void update_segments_var(int num_seg, Segment3D *segs, NodeType *nodeTypes, TimeLayer3D_GPU temp, TimeLayer3D_GPU next,  FTYPE *x, int  max_n_max_n, int id_shift = 0)
	{
		// fetch current segment info
		int id = id_shift + blockIdx.x * blockDim.x + threadIdx.x;
		if( id >= num_seg + id_shift) return;
		Segment3D &seg = segs[id];
		//NodeType &nodes = nodesBounds[id];

		switch (dir)
		{
		case X:
			if (seg.skipX)
				return;
		}
		update_segment<dir, var>(x, seg, nodeTypes, temp, next, id, num_seg, max_n_max_n);
	}

	template<int dir>
	void update_segments(dim3 grid, dim3 block, cudaStream_t &stream, int num_seg, Segment3D *segs, 
		                   NodeType *nodeTypes, TimeLayer3D_GPU temp, TimeLayer3D_GPU next,  FTYPE *d_x, int  max_n, int dimX, int id_shift = 0)
	{
		int max_n_max_n;
		int varOffset = (dimX + 2) * max_n * max_n * MAX_SEGS_PER_ROW;
		switch (dir)
		{
		case X:
			max_n_max_n = max_n * max_n;
			break;
		default:
			max_n_max_n = dimX * max_n; 
		}
		update_segments_var<dir, type_U><<<grid, block, 0, stream>>>(num_seg, segs, nodeTypes, temp, next, d_x + type_U*varOffset,  max_n_max_n, id_shift);
		update_segments_var<dir, type_V><<<grid, block, 0, stream>>>(num_seg, segs, nodeTypes, temp, next, d_x + type_V*varOffset,  max_n_max_n, id_shift);
		update_segments_var<dir, type_W><<<grid, block, 0, stream>>>(num_seg, segs, nodeTypes, temp, next, d_x + type_W*varOffset,  max_n_max_n, id_shift);
		update_segments_var<dir, type_T><<<grid, block, 0, stream>>>(num_seg, segs, nodeTypes, temp, next, d_x + type_T*varOffset,  max_n_max_n, id_shift);
	}

	template<int dir>
	void update_segments(dim3 grid, dim3 block, int num_seg, Segment3D *segs, 
		                   NodeType *nodeTypes, TimeLayer3D_GPU temp, TimeLayer3D_GPU next,  FTYPE *d_x, int  max_n, int dimX, int id_shift = 0)
	{
		cudaStream_t stream = 0;
		update_segments<dir>( grid, block, stream, num_seg, segs, nodeTypes, temp, next,  d_x, max_n, dimX, id_shift);
	}

	template<int dir, int swipe>
	void solve_segments(dim3 grid, dim3 block, cudaStream_t &stream, FluidParamsGPU **p, int num_seg, Segment3D *segs, NodesBoundary3D *nodesBounds, NodeType *nodeTypes, 
			                  TimeLayer3D_GPU &cur, TimeLayer3D_GPU &temp, TimeLayer3D_GPU &next, FTYPE *d_c, FTYPE *d_x, int max_n, int dimX, int id_shift = 0, int nBlockZ = 1)
	{
		int varOffset = ((max_n * max_n) % nBlockZ == 0) ? (max_n * max_n)/nBlockZ: (max_n * max_n)/nBlockZ + 1;
		varOffset *= (dimX + 2) * MAX_SEGS_PER_ROW;
		int max_n_max_n;
		switch (dir)
		{
		case X:
			max_n_max_n = max_n * max_n;
			break;
		case Y:
			max_n_max_n = dimX * max_n;
			max_n_max_n = ( max_n_max_n % nBlockZ == 0)? max_n_max_n / nBlockZ : max_n_max_n / nBlockZ + 1;
			break;
		default:
			max_n_max_n = dimX * max_n; 
		}

		solve_segments_var<dir, swipe, type_U><<<grid, block, 0, stream>>>( *p[type_U], num_seg, segs, nodesBounds, nodeTypes, cur, temp, next, d_c + type_U*varOffset, d_x + type_U*varOffset, max_n_max_n, dimX, id_shift );
		bool ok = (cudaSuccess == cudaGetLastError());
		solve_segments_var<dir, swipe, type_V><<<grid, block, 0, stream>>>( *p[type_V], num_seg, segs, nodesBounds, nodeTypes, cur, temp, next, d_c + type_V*varOffset, d_x + type_V*varOffset, max_n_max_n, dimX, id_shift );
		ok = (cudaSuccess == cudaGetLastError());
		solve_segments_var<dir, swipe, type_W><<<grid, block, 0, stream>>>( *p[type_W], num_seg, segs, nodesBounds, nodeTypes, cur, temp, next, d_c + type_W*varOffset, d_x + type_W*varOffset, max_n_max_n, dimX, id_shift );
		ok = (cudaSuccess == cudaGetLastError());
		solve_segments_var<dir, swipe, type_T><<<grid, block, 0, stream>>>( *p[type_T], num_seg, segs, nodesBounds, nodeTypes, cur, temp, next, d_c + type_T*varOffset, d_x + type_T*varOffset, max_n_max_n, dimX, id_shift );
		ok = (cudaSuccess == cudaGetLastError());
	}

	template<int dir, int swipe>
    void solve_segments(dim3 grid, dim3 block, FluidParamsGPU **p, int num_seg, Segment3D *segs, NodesBoundary3D *nodesBounds, NodeType *nodeTypes, 
			                  TimeLayer3D_GPU &cur, TimeLayer3D_GPU &temp, TimeLayer3D_GPU &next, FTYPE *d_c, FTYPE *d_x, int max_n, int dimX)
		{
			cudaStream_t stream = 0;
			solve_segments<dir, swipe>(grid, block, stream, p, num_seg, segs, nodesBounds, nodeTypes, cur, temp, next, d_c, d_x, max_n, dimX);
		}

	template<DirType dir>
	void LaunchSolveSegments_dir( FluidParamsGPU **p, int *num_seg, Segment3D **segs, NodesBoundary3D **nodesBounds, NodeType **nodeTypes, TimeLayer3D_GPU &cur, TimeLayer3D_GPU &temp, TimeLayer3D_GPU &next,
									  FTYPE **d_c, FTYPE **d_x )
// Y or Z direction only (and X if nGPUs = 1)
	{
		GPUplan *pGPUplan = GPUplan::Instance();
		PARAplan* pplan = PARAplan::Instance();

		int max_n = max( max( cur.dimx, cur.dimy ), cur.dimz );

		dim3 block(SOLVER_BLOCK_DIM);

		for (int i = 0; i < pGPUplan->size(); i++)
		{
			pGPUplan->setDevice(i);
			cur.SetDevice(i); temp.SetDevice(i); next.SetDevice(i);

			int dimX = pGPUplan->node(i)->getLength1D();
			dim3 grid((num_seg[i] + block.x - 1)/block.x);

			//cudaFuncSetCacheConfig(solve_segments<dir, var, ALL>, cudaFuncCachePreferL1);	 // will it create a queue?
            solve_segments<dir, ALL>(grid, block, p, num_seg[i], segs[i], nodesBounds[i], nodeTypes[i], cur, temp, next, d_c[i], d_x[i], max_n, dimX );
		}
		pGPUplan->deviceSynchronize();
	}

	//template<VarType var>
	void LaunchSolveSegments_X( FluidParamsGPU **p, int *num_seg, Segment3D **segs, NodesBoundary3D **nodesBounds, NodeType **nodeTypes, TimeLayer3D_GPU &cur, TimeLayer3D_GPU &temp, TimeLayer3D_GPU &next,
									  FTYPE **d_c, FTYPE **d_x, int numSegs, FTYPE *mpi_buf = NULL )
	/*
		X direction only
	*/
	{
		int typeID;
		int varOffset;
		GPUplan *pGPUplan = GPUplan::Instance();
		PARAplan* pplan = PARAplan::Instance();
		int irank =  pplan->rank();
		int size = pplan->size();        
		
		if (pGPUplan->size() == 1 && pplan->size() == 1)
		{
			LaunchSolveSegments_dir<X>( p, num_seg, segs, nodesBounds, nodeTypes, cur, temp, next, d_c, d_x );
			return;
		}
		
		int max_n = max( max( cur.dimx, cur.dimy ), cur.dimz );

		int haloSize = max_n * max_n * MAX_SEGS_PER_ROW;
		int comSize = numSegs;
		int dimX;

		dim3 block(SOLVER_BLOCK_DIM);
#ifdef __PARA
		if (pplan->size() > 1)
		{	
			pGPUplan->setDevice(0);
			dimX = pGPUplan->node(0)->getLength1D();
			varOffset = (dimX + 2) * max_n * max_n * MAX_SEGS_PER_ROW;
			for (typeID = 0; typeID < SOLVER_VAR_NUM; typeID++)
			{
				paraDevRecv<FTYPE, FORWARD>(d_c[0] + typeID*varOffset, mpi_buf, comSize, 666 + typeID);
				paraDevRecv<FTYPE, FORWARD>(d_x[0] + typeID*varOffset, mpi_buf + comSize, comSize, 666 + SOLVER_VAR_NUM + typeID);
			}
		}
#endif
		for (int i = 0; i < pGPUplan->size(); i++)
		{
			pGPUplan->setDevice(i);
			cur.SetDevice(i); temp.SetDevice(i); next.SetDevice(i);

			dimX = pGPUplan->node(i)->getLength1D();
			dim3 grid((num_seg[i] + block.x - 1)/block.x);
			varOffset = (dimX + 2) * max_n * max_n * MAX_SEGS_PER_ROW;
			//cudaFuncSetCacheConfig(solve_segments<dir, var, ALL>, cudaFuncCachePreferL1);
            solve_segments<X, FORWARD>(grid, block, p, num_seg[i], segs[i], nodesBounds[i], nodeTypes[i], cur, temp, next, d_c[i], d_x[i], max_n, dimX );

			for (typeID = 0; typeID < SOLVER_VAR_NUM; typeID++)
				if (i < pGPUplan->size() - 1) //  send to node n+1
				{	
					int varOffset_p1 = (pGPUplan->node(i+1)->getLength1D() + 2) * max_n * max_n * MAX_SEGS_PER_ROW;
					haloMemcpyPeer<FTYPE, FORWARD>( d_c, i, haloSize, dimX * haloSize + typeID*(varOffset - varOffset_p1), typeID*varOffset_p1, comSize);					
					haloMemcpyPeer<FTYPE, FORWARD>( d_x, i, haloSize, dimX * haloSize + typeID*(varOffset - varOffset_p1), typeID*varOffset_p1, comSize);
				}			
		}
#ifdef __PARA
		if (pplan->size() > 0)
		{
			pGPUplan->setDevice(pGPUplan->size()-1);
			dimX = pGPUplan->node(pGPUplan->size()-1)->getLength1D();
			varOffset = (dimX + 2) * max_n * max_n * MAX_SEGS_PER_ROW;
			for (typeID = 0; typeID < SOLVER_VAR_NUM; typeID++)
			{
				paraDevSend<FTYPE, FORWARD>(d_c[pGPUplan->size()-1] + typeID * varOffset + haloSize + dimX * haloSize - haloSize, mpi_buf, comSize, 666 + typeID);
				paraDevSend<FTYPE, FORWARD>(d_x[pGPUplan->size()-1]  + typeID * varOffset + haloSize + dimX * haloSize - haloSize, mpi_buf + comSize, comSize, 666 + SOLVER_VAR_NUM + typeID);				
			}
			for (typeID = 0; typeID < SOLVER_VAR_NUM; typeID++)
				paraDevRecv<FTYPE, BACK>(d_x[pGPUplan->size()-1] + typeID * varOffset + haloSize +  dimX * haloSize, mpi_buf, comSize, 666 + 2 * SOLVER_VAR_NUM + typeID);
		}
#endif
		for (int i = pGPUplan->size() - 1; i >= 0; i--)
		{
			pGPUplan->setDevice(i);
			cur.SetDevice(i); temp.SetDevice(i); next.SetDevice(i);

			dimX = pGPUplan->node(i)->getLength1D();
			dim3 grid((num_seg[i] + block.x - 1)/block.x);
			varOffset = (dimX + 2) * max_n * max_n * MAX_SEGS_PER_ROW;
			solve_segments<X, BACK>(grid, block, p, num_seg[i], segs[i], nodesBounds[i], nodeTypes[i], cur, temp, next, d_c[i], d_x[i], max_n, dimX );
			for (typeID = 0; typeID < SOLVER_VAR_NUM; typeID++)
				if (i > 0)
				{
					int dimXm1 = pGPUplan->node(i-1)->getLength1D();
					int varOffset_m1 = (dimXm1 + 2) * max_n * max_n * MAX_SEGS_PER_ROW;					
					haloMemcpyPeer<FTYPE, BACK>(d_x, i, haloSize, dimXm1*haloSize + typeID * (varOffset_m1 - varOffset), typeID*varOffset, comSize);
				}
		}
#ifdef __PARA		
		if (pplan->size() > 0)
		{			
			pGPUplan->setDevice(0);
			dimX = pGPUplan->node(0)->getLength1D();
			varOffset = (dimX + 2) * max_n * max_n * MAX_SEGS_PER_ROW;
			for (typeID = 0; typeID < SOLVER_VAR_NUM; typeID++)
				paraDevSend<FTYPE, BACK>(d_x[0] + typeID * varOffset + haloSize, mpi_buf, comSize, 666 + 2 * SOLVER_VAR_NUM + typeID);
		}
#endif
		for (int i = 0; i < pGPUplan->size(); i++)
		{
			pGPUplan->setDevice(i);
			cur.SetDevice(i); temp.SetDevice(i); next.SetDevice(i);
			dim3 grid((num_seg[i] + block.x - 1)/block.x);
			//for (int typeID = 0; typeID < SOLVER_VAR_NUM; typeID++)
			dimX = pGPUplan->node(i)->getLength1D();
			varOffset = (dimX + 2) * max_n * max_n * MAX_SEGS_PER_ROW;
			update_segments<X>(grid, block, num_seg[i], segs[i], nodeTypes[i], temp, next, d_x[i], max_n, dimX);
			//update_segments_var<X, type_U><<<grid, block>>>(num_seg[i], segs[i], nodeTypes[i], temp, next, d_x[i] + type_U*varOffset,  max_n_max_n);
			//update_segments_var<X, type_V><<<grid, block>>>(num_seg[i], segs[i], nodeTypes[i], temp, next, d_x[i] + type_V*varOffset,  max_n_max_n);
			//update_segments_var<X, type_W><<<grid, block>>>(num_seg[i], segs[i], nodeTypes[i], temp, next, d_x[i] + type_W*varOffset,  max_n_max_n);
			//update_segments_var<X, type_T><<<grid, block>>>(num_seg[i], segs[i], nodeTypes[i], temp, next, d_x[i] + type_T*varOffset,  max_n_max_n);
			//update_segments<X, var><<<grid, block>>>( num_seg[i], segs[i], nodeTypes[i], temp, next, d_x[i],  max_n_max_n);
		}
		pGPUplan->deviceSynchronize();
	}

	void LaunchSolveSegments_XY( FluidParamsGPU **px, FluidParamsGPU **py, int **num_segXBlZ,  int **num_segYBlZ, int **comuNumSegsXBlZ, int **comuNumSegsYBlZ, Segment3D **segsX, Segment3D **segsY, 
		                                   int num_local, int nblock, NodesBoundary3D **nodesBoundsX, NodesBoundary3D **nodesBoundsY, NodeType **nodeTypes, 
																			 TimeLayer3D_GPU &cur, TimeLayer3D_GPU &temp, TimeLayer3D_GPU &half, TimeLayer3D_GPU &next,
																			 FTYPE **d_c, FTYPE **d_cY, FTYPE **d_x, FTYPE **d_xY )
	{ 
		PARAplan* pplan = PARAplan::Instance();
#ifdef __PARA
		if (pplan->size() > 1)
			throw std::logic_error("Not Implemented");
#endif

#if !INTERNAL_MERGE_ENABLE
	throw std::logic_error("Not Implemented");
#endif

		int typeID, varOffset;
		int max_n = max( max( cur.dimx, cur.dimy ), cur.dimz );
		int max_n_max_n_X = max_n * max_n;

		// Need to sync halos at each iteration!!!

		GPUplan *pGPUplan = GPUplan::Instance();
		dim3 block(SOLVER_BLOCK_DIM);
		int nnodes =  pGPUplan->size();
		int haloSize, dimX;

		int i, it, iblock, iblock_it, iblock_back_it;

		for (iblock = 0 ; iblock < nblock + (2 * (nnodes - 1) + 1) * num_local - 1; iblock ++)
		{
			// Solve Y segment of iblock
			if (iblock < nblock)
				for (int it = 0; it < num_local; it++)
				{
					for ( i = 0; i < nnodes; i++)		
					{
						cudaSetDevice(i);
						cur.SetDevice(i); temp.SetDevice(i); half.SetDevice(i); next.SetDevice(i);
						int dimX = pGPUplan->node(i)->getLength1D();
						//max_n_max_n_Y = max_n * pGPUplan->node(i)->getLength1D();
						//max_n_max_n_Y = ( max_n_max_n_Y % nblock == 0)? max_n_max_n_Y / nblock : max_n_max_n_Y / nblock + 1;
						dim3 grid((num_segYBlZ[i][iblock] + block.x - 1)/block.x);
						//segsY[i] + comuNumSegsYBlZ[i][iblock]
						//solve_segments<Y, ALL>(grid, block,  pGPUplan->node(i)->stream, py, num_segYBlZ[i][iblock] + comuNumSegsYBlZ[i][iblock], segsY[i], 
						//												 nodesBoundsY[i], nodeTypes[i], cur, temp, half, d_cY[i],  d_xY[i], max_n, dimX, 0, nblock); //???
						solve_segments<Y, ALL>(grid, block,  pGPUplan->node(i)->stream, py, num_segYBlZ[i][iblock], segsY[i], 
																		 nodesBoundsY[i], nodeTypes[i], cur, temp, half, d_cY[i],  d_xY[i], max_n, dimX, comuNumSegsYBlZ[i][iblock], nblock); //???
					}
					for ( i = 0; i < nnodes; i++)
					{
						pGPUplan->setDevice(i);
						gpuSafeCall(cudaStreamSynchronize(pGPUplan->node(i)->stream), "multiDevHaloSync(): cudaStreamSynchronize", i, __FILE__, __LINE__);
					}
					//sync:					
					haloSize = cur.dimy * cur.dimz;
					multiDevHaloSync<FTYPE>(temp.dd_u, haloSize, false, iblock*haloSize / nblock, haloSize / nblock);
					//throw std::logic_error("Testing");
					multiDevHaloSync<FTYPE>(temp.dd_v, haloSize, false, iblock*haloSize / nblock, haloSize / nblock);
					multiDevHaloSync<FTYPE>(temp.dd_w, haloSize, false, iblock*haloSize / nblock, haloSize / nblock);
					multiDevHaloSync<FTYPE>(temp.dd_T, haloSize, true, iblock*haloSize / nblock, haloSize / nblock);						
				}

				/*for ( i = 0; i < nnodes; i++)
				{
					pGPUplan->setDevice(i);
					gpuSafeCall(cudaStreamSynchronize(pGPUplan->node(i)->stream), "multiDevHaloSync(): cudaStreamSynchronize", i, __FILE__, __LINE__);
					gpuSafeCall(cudaStreamSynchronize(pGPUplan->node(i)->stream2), "multiDevHaloSync(): cudaStreamSynchronize", i, __FILE__, __LINE__);
				}*/
        
			// Sync complete iterations:
			//for (it = 1; it < num_local; it++)
			//{
			//	//iblock_it = (nnodes > 2)? iblock - 2 * it * (nnodes-1): iblock - 2 * it; 
			//	iblock_it = iblock - 2 * it * (nnodes-1);
			//	if (iblock_it < nblock && iblock_it >= 0)
			//	{
			//		haloSize = cur.dimy * cur.dimz;			
			//		multiDevHaloSync<FTYPE>(temp.dd_u, haloSize, false, iblock_it*haloSize / nblock, haloSize / nblock);
			//		multiDevHaloSync<FTYPE>(temp.dd_v, haloSize, false, iblock_it*haloSize / nblock, haloSize / nblock);
			//		multiDevHaloSync<FTYPE>(temp.dd_w, haloSize, false, iblock_it*haloSize / nblock, haloSize / nblock);
			//		multiDevHaloSync<FTYPE>(temp.dd_T, haloSize, true, iblock_it*haloSize / nblock, haloSize / nblock);
			//	}
			//}
			//Wait for incoming X data from node i-1
			for (i = 1; i < nnodes; i++)
			{
				cudaSetDevice(i);
				for (it = 0; it < num_local; it++)
				{
					iblock_it = iblock -  it * (2 * (nnodes-1) + 1) - i; 
					if (iblock_it >= 0 && i > 0 && iblock_it < nblock)	
						gpuSafeCall( cudaStreamWaitEvent(pGPUplan->node(i)->stream, pGPUplan->node(i-1)->event[2 * iblock_it], 0), "LaunchSolveSegments_XY: cudaStreamWaitEvent" ); //event[0]
				}
			}

			// Solve X segments forward
			for (it = 0; it < num_local; it++)			
				for (i = 0; i < nnodes; i++)
				{ 
					iblock_it = iblock - it * (2 * (nnodes-1) + 1) - i;
					if (iblock_it < nblock && iblock_it >= 0)
					{
						cudaSetDevice(i);
						cur.SetDevice(i); temp.SetDevice(i); half.SetDevice(i); next.SetDevice(i);
						dimX = pGPUplan->node(i)->getLength1D();
						//int num_elems = dimX * max_n_max_n_X * MAX_SEGS_PER_ROW;

						//cudaStreamSynchronize(pGPUplan->node(i)->stream); // sync with Y calculation
						dim3 grid(( num_segXBlZ[i][iblock_it] + block.x - 1)/block.x);
						//solve_segments<X, FORWARD><<<grid, block, 0, pGPUplan->node(i)->stream>>>( num_segXBlZ[i][iblock_it], segsX[i], 
						//															nodes[i], half, temp, next, dd_a[i],  dd_x[i], max_n_max_n_X, dimX, comuNumSegsXBlZ[i][iblock_it]);
						solve_segments<X, FORWARD>(grid, block,  pGPUplan->node(i)->stream, px, num_segXBlZ[i][iblock_it], segsX[i], 
																		nodesBoundsX[i], nodeTypes[i], half, temp, next, d_c[i],  d_x[i], max_n, dimX, comuNumSegsXBlZ[i][iblock_it]);						
					}
				}
			
			//Send data to node i+1
			for (i = 0; i < nnodes - 1; i++)
			{
				haloSize = max_n_max_n_X * MAX_SEGS_PER_ROW;
				dimX = pGPUplan->node(i)->getLength1D();
				varOffset = (dimX + 2) * max_n * max_n * MAX_SEGS_PER_ROW;
				int varOffset_p1 = (pGPUplan->node(i+1)->getLength1D() + 2) * max_n * max_n * MAX_SEGS_PER_ROW;
				for (int it = 0; it < num_local; it++)
				{
					iblock_it = iblock - it * (2 * (nnodes-1) + 1) - i; 
					if (iblock_it < nblock && iblock_it >= 0)
					{						
						for (typeID = 0; typeID < SOLVER_VAR_NUM; typeID++)
						{
							//haloMemcpyPeerAsync<FTYPE, FORWARD>(d_c, i, haloSize, dimX * haloSize, pGPUplan->node(i)->stream, comuNumSegsXBlZ[i][iblock_it], num_segXBlZ[i][iblock_it]);							
							haloMemcpyPeerAsync<FTYPE, FORWARD>( d_c, i, haloSize, dimX * haloSize + typeID*(varOffset - varOffset_p1), pGPUplan->node(i)->stream, typeID*varOffset_p1 + comuNumSegsXBlZ[i][iblock_it], num_segXBlZ[i][iblock_it]);					
							haloMemcpyPeerAsync<FTYPE, FORWARD>( d_x, i, haloSize, dimX * haloSize + typeID*(varOffset - varOffset_p1), pGPUplan->node(i)->stream, typeID*varOffset_p1 + comuNumSegsXBlZ[i][iblock_it], num_segXBlZ[i][iblock_it]);
						}
						cudaSetDevice(i);
						gpuSafeCall( cudaEventRecord( pGPUplan->node(i)->event[2 * iblock_it],  pGPUplan->node(i)->stream), "");
					}
				}
			}
			
			cudaStreamSynchronize(pGPUplan->node(nnodes-1)->stream); // sync for immediate back computation on last device
			
			//Wait for incoming X data from node i+1
			for (i = 0; i < nnodes  - 1; i++)
			{
				cudaSetDevice(i);
				for (it = 0; it < num_local; it++)
				{
					iblock_back_it = iblock - 2 * (nnodes-1) - it * (2 * (nnodes-1) + 1) + i;
					if (iblock_back_it < nblock && iblock_back_it >= 0)
						gpuSafeCall( cudaStreamWaitEvent(pGPUplan->node(i)->stream2, pGPUplan->node(i+1)->event[2 * iblock_back_it + 1], 0), ""); // event[1]
				}
			}
			// Solve X segments backwards
			for (it = 0; it < num_local; it++)		
			{
				for (i = 0; i < nnodes; i++)
				{
					cudaSetDevice(i);
					cur.SetDevice(i); temp.SetDevice(i); half.SetDevice(i); next.SetDevice(i);
					dimX = pGPUplan->node(i)->getLength1D();
					varOffset = (dimX + 2) * max_n * max_n * MAX_SEGS_PER_ROW;				
					iblock_back_it = iblock - 2 * (nnodes-1) - it * (2 * (nnodes-1) + 1) + i;
					if (iblock_back_it < nblock && iblock_back_it >= 0)
					{
						dim3 grid(( num_segXBlZ[i][iblock_back_it] + block.x - 1)/block.x);
						//solve_segments<X, BACK><<<grid, block, 0, pGPUplan->node(i)->stream2>>>( num_segXBlZ[i][iblock_back_it], segsX[i], 
  					 //																						nodes[i], half, temp, next, dd_a[i],  dd_x[i], max_n_max_n_X, dimX, comuNumSegsXBlZ[i][iblock_back_it]);
						solve_segments<X, BACK>(grid, block,  pGPUplan->node(i)->stream2, px, num_segXBlZ[i][iblock_back_it], segsX[i], 
																	 nodesBoundsX[i], nodeTypes[i], half, temp, next, d_c[i],  d_x[i], max_n, dimX, comuNumSegsXBlZ[i][iblock_back_it]);
						if (num_local > 1)
							update_segments<X>(grid, block, pGPUplan->node(i)->stream2, num_segXBlZ[i][iblock_back_it], segsX[i], nodeTypes[i], temp, next, d_x[i], max_n, dimX, comuNumSegsXBlZ[i][iblock_back_it]);
					}
						//************************
				}
				// put nonlinear sync here:
				if (num_local > 1) {
					iblock_back_it = iblock - 2 * (nnodes-1) - it * (2 * (nnodes-1) + 1);
					if (iblock_back_it < nblock && iblock_back_it >= 0)
					{
						haloSize = cur.dimy * cur.dimz;						
						multiDevHaloSync2<FTYPE>(temp.dd_u, haloSize, false, iblock_back_it*haloSize / nblock, haloSize / nblock);
						multiDevHaloSync2<FTYPE>(temp.dd_v, haloSize, false, iblock_back_it*haloSize / nblock, haloSize / nblock);
						multiDevHaloSync2<FTYPE>(temp.dd_w, haloSize, false, iblock_back_it*haloSize / nblock, haloSize / nblock);
						multiDevHaloSync2<FTYPE>(temp.dd_T, haloSize, true, iblock_back_it*haloSize / nblock, haloSize / nblock);			
					}	
				}
			}

			//Send data to node i-1
			for (i = 1; i < nnodes; i++)
			{
				haloSize = max_n_max_n_X * MAX_SEGS_PER_ROW;
				dimX = pGPUplan->node(i)->getLength1D();
				int dimXm1 = pGPUplan->node(i-1)->getLength1D();
				int varOffset_m1 = (dimXm1 + 2) * max_n * max_n * MAX_SEGS_PER_ROW;
				for (int it = 0; it < num_local; it++)
				{
					iblock_back_it = iblock - 2 * (nnodes-1) - it * (2 * (nnodes-1) + 1) + i;
					if (iblock_back_it < nblock && iblock_back_it >= 0)
					{
						for (typeID = 0; typeID < SOLVER_VAR_NUM; typeID++)
						{
							//haloMemcpyPeerAsync<FTYPE, BACK>(d_x, i, haloSize, pGPUplan->node(i-1)->getLength1D() * haloSize, pGPUplan->node(i)->stream2, comuNumSegsXBlZ[i-1][iblock_back_it], num_segXBlZ[i][iblock_back_it]);														
							haloMemcpyPeerAsync<FTYPE, BACK>(d_x, i, haloSize, dimXm1*haloSize + typeID * (varOffset_m1 - varOffset), pGPUplan->node(i)->stream2, typeID*varOffset + comuNumSegsXBlZ[i-1][iblock_back_it], num_segXBlZ[i][iblock_back_it]);
						}
						cudaSetDevice(i);
						gpuSafeCall( cudaEventRecord( pGPUplan->node(i)->event[2 * iblock_back_it + 1],  pGPUplan->node(i)->stream2), "" );
					}
				}
			}
			// put nonlinear sync here:
			if (num_local == 1) {
				for (i = 0; i < nnodes; i++)
				{
					cudaSetDevice(i);
					cur.SetDevice(i); temp.SetDevice(i); half.SetDevice(i); next.SetDevice(i);
					dimX = pGPUplan->node(i)->getLength1D();
					varOffset = (dimX + 2) * max_n * max_n * MAX_SEGS_PER_ROW;				
					iblock_back_it = iblock - 2 * (nnodes-1) + i;
					if (iblock_back_it < nblock && iblock_back_it >= 0)
					{
						dim3 grid(( num_segXBlZ[i][iblock_back_it] + block.x - 1)/block.x);
						update_segments<X>(grid, block, pGPUplan->node(i)->stream2, num_segXBlZ[i][iblock_back_it], segsX[i], nodeTypes[i], temp, next, d_x[i], max_n, dimX, comuNumSegsXBlZ[i][iblock_back_it]);						
					}
				}					
				iblock_back_it = iblock - 2 * (nnodes-1);
				if (iblock_back_it < nblock && iblock_back_it >= 0)
				{
					haloSize = cur.dimy * cur.dimz;						
					multiDevHaloSync2<FTYPE>(temp.dd_u, haloSize, false, iblock_back_it*haloSize / nblock, haloSize / nblock);
					multiDevHaloSync2<FTYPE>(temp.dd_v, haloSize, false, iblock_back_it*haloSize / nblock, haloSize / nblock);
					multiDevHaloSync2<FTYPE>(temp.dd_w, haloSize, false, iblock_back_it*haloSize / nblock, haloSize / nblock);
					multiDevHaloSync2<FTYPE>(temp.dd_T, haloSize, true, iblock_back_it*haloSize / nblock, haloSize / nblock);			
				}				
			}
		}
		//throw std::logic_error("Testing...");
		pGPUplan->deviceSynchronize();
	}
		
	//void SolveSegments_GPU( int* num_seg, Segment3D **segs, DirType dir, Node **nodes, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *next,
	//						FTYPE **dd_a, FTYPE **dd_x )
	//{
	//	TimeLayer3D_GPU d_cur( cur );
	//	TimeLayer3D_GPU d_temp( temp );
	//	TimeLayer3D_GPU d_next( next );

	//	switch( dir )
	//	{
	//	case X: LaunchSolveSegments_X( num_seg, segs, nodes, d_cur, d_temp, d_next, dd_a, dd_x ); break;
	//	case Y: LaunchSolveSegments_dir<Y>( num_seg, segs, nodes, d_cur, d_temp, d_next, dd_a, dd_x ); break;
	//	case Z: LaunchSolveSegments_dir<Z>( num_seg, segs, nodes, d_cur, d_temp, d_next, dd_a, dd_x ); break;
	//	}
	//}

	//template<DirType dir>
	//void LaunchSolveSegments_dir( FluidParamsGPU p, int *num_seg, Segment3D **segs, VarType var, NodesBoundary3D **nodesBounds, NodeType **nodeTypes, TimeLayer3D_GPU &cur, TimeLayer3D_GPU &temp, TimeLayer3D_GPU &next,
	//							  FTYPE **d_c, FTYPE **d_x )
	//{
	//	switch( var )
	//	{
	//	case type_U: LaunchSolveSegments_dir_var<dir, type_U>( p, num_seg, segs, nodesBounds, nodeTypes, cur, temp, next, d_c, d_x ); break;
	//	case type_V: LaunchSolveSegments_dir_var<dir, type_V>( p, num_seg, segs, nodesBounds, nodeTypes, cur, temp, next, d_c, d_x ); break;
	//	case type_W: LaunchSolveSegments_dir_var<dir, type_W>( p, num_seg, segs, nodesBounds, nodeTypes, cur, temp, next, d_c, d_x ); break;
	//	case type_T: LaunchSolveSegments_dir_var<dir, type_T>( p, num_seg, segs, nodesBounds, nodeTypes, cur, temp, next, d_c, d_x ); break;
	//	}
	//}

	//void LaunchSolveSegments_X( FluidParamsGPU p, int *num_seg, Segment3D **segs, VarType var, NodesBoundary3D **nodesBounds, NodeType **nodeTypes, TimeLayer3D_GPU &cur, TimeLayer3D_GPU &temp, TimeLayer3D_GPU &next,
	//							  FTYPE **d_c, FTYPE **d_x, int numSegs, FTYPE *mpi_buf = NULL )
	//{
	//	switch( var )
	//	{
	//	case type_U: LaunchSolveSegments_X_var<type_U>( p, num_seg, segs, nodesBounds, nodeTypes, cur, temp, next, d_c, d_x, numSegs, mpi_buf ); break;
	//	case type_V: LaunchSolveSegments_X_var<type_V>( p, num_seg, segs, nodesBounds, nodeTypes, cur, temp, next, d_c, d_x, numSegs, mpi_buf ); break;
	//	case type_W: LaunchSolveSegments_X_var<type_W>( p, num_seg, segs, nodesBounds, nodeTypes, cur, temp, next, d_c, d_x, numSegs, mpi_buf ); break;
	//	case type_T: LaunchSolveSegments_X_var<type_T>( p, num_seg, segs, nodesBounds, nodeTypes, cur, temp, next, d_c, d_x, numSegs, mpi_buf ); break;
	//	}
	//}

	void SolveSegments_GPU( FTYPE dt, FluidParams params, int *num_seg, Segment3D **segs, DirType dir, NodesBoundary3D **nodesBounds,  NodeType **nodeTypes, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *next,
							FTYPE **d_c, FTYPE **d_x,  int numSegs, FTYPE *mpi_buf )
	{
		TimeLayer3D_GPU d_cur( cur );
		TimeLayer3D_GPU d_temp( temp );
		TimeLayer3D_GPU d_next( next );

		FluidParamsGPU *p[SOLVER_VAR_NUM];

		for (int ip = 0; ip < SOLVER_VAR_NUM; ip++)
		{
			p[ip] = new FluidParamsGPU( VarType(ip), dir, dt, cur->dx, cur->dy, cur->dz, params );
		}
		//FluidParamsGPU pp( var, dir, dt, cur->dx, cur->dy, cur->dz, params );

		switch( dir )
		{
		case X: LaunchSolveSegments_X( p, num_seg, segs, nodesBounds, nodeTypes, d_cur, d_temp, d_next, d_c, d_x, numSegs, mpi_buf ); break;
		case Y: LaunchSolveSegments_dir<Y>( p, num_seg, segs, nodesBounds, nodeTypes, d_cur, d_temp, d_next, d_c, d_x ); break;
		case Z: LaunchSolveSegments_dir<Z>( p, num_seg, segs, nodesBounds, nodeTypes, d_cur, d_temp, d_next, d_c, d_x ); break;
		case Z_as_Y: LaunchSolveSegments_dir<Z_as_Y>( p, num_seg, segs, nodesBounds, nodeTypes, d_cur, d_temp, d_next, d_c, d_x ); break;
		}

		for (int ip = 0; ip < SOLVER_VAR_NUM; ip ++)
			delete p[ip];
	}

	void SolveSegments_XY_GPU( FTYPE dt, FluidParams params, int **num_segXBlZ,  int **num_segYBlZ, int **comuNumSegsXBlZ, int **comuNumSegsYBlZ, Segment3D **segsX, Segment3D **segsY, 
		                                   int num_local, int nblock, NodesBoundary3D **nodesBoundsX, NodesBoundary3D **nodesBoundsY, NodeType **nodeTypes, 
																			 TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *half, TimeLayer3D *next,
																			 FTYPE **d_c, FTYPE **d_cY, FTYPE **d_x, FTYPE **d_xY )
	{
		TimeLayer3D_GPU d_cur( cur );
		TimeLayer3D_GPU d_temp( temp );
		TimeLayer3D_GPU d_next( next );
		TimeLayer3D_GPU d_half( half );

		FluidParamsGPU *px[SOLVER_VAR_NUM];
		FluidParamsGPU *py[SOLVER_VAR_NUM];

		for (int ip = 0; ip < SOLVER_VAR_NUM; ip++)
		{
			px[ip] = new FluidParamsGPU( VarType(ip), X, dt, cur->dx, cur->dy, cur->dz, params );
			py[ip] = new FluidParamsGPU( VarType(ip), Y, dt, cur->dx, cur->dy, cur->dz, params );
		}

		LaunchSolveSegments_XY(px, py, num_segXBlZ, num_segYBlZ, comuNumSegsXBlZ, comuNumSegsYBlZ, segsX, segsY, num_local, 
			                      nblock, nodesBoundsX, nodesBoundsY, nodeTypes, d_cur, d_temp, d_half, d_next, d_c, d_cY, d_x, d_xY);		

		for (int ip = 0; ip < SOLVER_VAR_NUM; ip ++)
		{
			delete py[ip];
			delete px[ip];
		}
	}
}
