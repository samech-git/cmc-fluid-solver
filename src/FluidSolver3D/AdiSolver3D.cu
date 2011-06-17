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

#define SOLVER_BLOCK_DIM	256

#define SEG_BLOCK_DIM_X		32
#define SEG_BLOCK_DIM_Y		8

namespace FluidSolver3D
{
	struct FluidParamsGPU
	{
		FTYPE vis_dx2;
		FTYPE dt, dx, dy, dz;
		FTYPE v_T, t_phi;

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
	#define get(a, elem_id)			a[id + (elem_id) * max_n * max_n * MAX_SEGS_PER_ROW]
#else
	// sequential layout - bad access pattern
	#define get(a, elem_id)			a[elem_id + id * max_n]
#endif

	template<int dir, int var>
	__device__ void apply_bc0(int i, int j, int k, int dimy, int dimz, FTYPE &b0, FTYPE &c0, FTYPE &d0, Node *nodes)
	{
		int id = i * dimy * dimz;
		if (dir != Z_as_Y) id += j * dimz + k;
			else id += j * dimy + k;

		if ((var == type_T && nodes[id].bc_temp == BC_FREE) ||
			(var != type_T && nodes[id].bc_vel == BC_FREE))
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
			case type_U: d0 = nodes[id].v.x; break;
			case type_V: d0 = nodes[id].v.y; break;
			case type_W: d0 = nodes[id].v.z; break;
			case type_T: d0 = nodes[id].T; break;
			}
		}
	}

	template<int dir, int var>
	__device__ void apply_bc1(int i, int j, int k, int dimy, int dimz, FTYPE &a1, FTYPE &b1, FTYPE &d1, Node *nodes)
	{
		int id = i * dimy * dimz;
		if (dir != Z_as_Y) id += j * dimz + k;
			else id += j * dimy + k;

		if ((var == type_T && nodes[id].bc_temp == BC_FREE) ||
			(var != type_T && nodes[id].bc_vel == BC_FREE))
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
			case type_U: d1 = nodes[id].v.x; break;
			case type_V: d1 = nodes[id].v.y; break;
			case type_W: d1 = nodes[id].v.z; break;
			case type_T: d1 = nodes[id].T; break;
			}
		}
	}	

	template<int dir, int var>
	__device__ void build_matrix(FluidParamsGPU params, int i, int j, int k, FTYPE *a, FTYPE *b, FTYPE *c, FTYPE *d, int n, TimeLayer3D_GPU &cur, TimeLayer3D_GPU &temp, int id, int num_seg, int max_n)
	{	
		for (int p = 1; p < n-1; p++)
		{
			switch (dir)
			{
			case X:		
				get(a,p) = - temp.elem(temp.u, i+p, j, k) / (2 * params.dx) - params.vis_dx2; 
				get(b,p) = 3 / params.dt  +  2 * params.vis_dx2; 
				get(c,p) = temp.elem(temp.u, i+p, j, k) / (2 * params.dx) - params.vis_dx2; 
				
				switch (var)	
				{
				case type_U: get(d,p) = cur.elem(cur.u, i+p, j, k) * 3 / params.dt - params.v_T * temp.d_x(temp.T, params.dx, i+p, j, k); break;
				case type_V: get(d,p) = cur.elem(cur.v, i+p, j, k) * 3 / params.dt; break;
				case type_W: get(d,p) = cur.elem(cur.w, i+p, j, k) * 3 / params.dt; break;
				case type_T: get(d,p) = cur.elem(cur.T, i+p, j, k) * 3 / params.dt + params.t_phi * temp.DissFuncX(params.dx, params.dy, params.dz, i+p, j, k); break;
				}	
				break;

			case Y:
				get(a,p) = - temp.elem(temp.v, i, j+p, k) / (2 * params.dy) - params.vis_dx2; 
				get(b,p) = 3 / params.dt  +  2 * params.vis_dx2; 
				get(c,p) = temp.elem(temp.v, i, j+p, k) / (2 * params.dy) - params.vis_dx2; 
				
				switch (var)	
				{
				case type_U: get(d,p) = cur.elem(cur.u, i, j+p, k) * 3 / params.dt; break;
				case type_V: get(d,p) = cur.elem(cur.v, i, j+p, k) * 3 / params.dt - params.v_T * temp.d_y(temp.T, params.dy, i, j+p, k); break;
				case type_W: get(d,p) = cur.elem(cur.w, i, j+p, k) * 3 / params.dt; break;
				case type_T: get(d,p) = cur.elem(cur.T, i, j+p, k) * 3 / params.dt + params.t_phi * temp.DissFuncY(params.dx, params.dy, params.dz, i, j+p, k); break;
				}
				break;

			case Z:
				get(a,p) = - temp.elem(temp.w, i, j, k+p) / (2 * params.dz) - params.vis_dx2; 
				get(b,p) = 3 / params.dt  +  2 * params.vis_dx2; 
				get(c,p) = temp.elem(temp.w, i, j, k+p) / (2 * params.dz) - params.vis_dx2; 
				
				switch (var)	
				{
				case type_U: get(d,p) = cur.elem(cur.u, i, j, k+p) * 3 / params.dt; break;
				case type_V: get(d,p) = cur.elem(cur.v, i, j, k+p) * 3 / params.dt; break;
				case type_W: get(d,p) = cur.elem(cur.w, i, j, k+p) * 3 / params.dt - params.v_T * temp.d_z(temp.T, params.dz, i, j, k+p); break;
				case type_T: get(d,p) = cur.elem(cur.T, i, j, k+p) * 3 / params.dt + params.t_phi * temp.DissFuncZ(params.dx, params.dy, params.dz, i, j, k+p); break;
				}
				break;

			case Z_as_Y:
				get(a,p) = - temp.elem(temp.w, i, k+p, j) / (2 * params.dz) - params.vis_dx2; 
				get(b,p) = 3 / params.dt  +  2 * params.vis_dx2; 
				get(c,p) = temp.elem(temp.w, i, k+p, j) / (2 * params.dz) - params.vis_dx2; 
				
				switch (var)	
				{
				case type_U: get(d,p) = cur.elem(cur.u, i, k+p, j) * 3 / params.dt; break;
				case type_V: get(d,p) = cur.elem(cur.v, i, k+p, j) * 3 / params.dt; break;
				case type_W: get(d,p) = cur.elem(cur.w, i, k+p, j) * 3 / params.dt - params.v_T * temp.d_y(temp.T, params.dz, i, k+p, j); break;
				case type_T: get(d,p) = cur.elem(cur.T, i, k+p, j) * 3 / params.dt + params.t_phi * temp.DissFuncZ_as_Y(params.dx, params.dz, params.dy, i, k+p, j); break;
				}
				break;
			}
		}
	}

	__device__ void solve_tridiagonal( FTYPE *a, FTYPE *b, FTYPE *c, FTYPE *d, FTYPE *x, int num, int id, int num_seg, int max_n )
	{
		get(c,num-1) = 0.0;
		
		get(c,0) = get(c,0) / get(b,0);
		get(d,0) = get(d,0) / get(b,0);

		for (int i = 1; i < num; i++)
		{
			get(c,i) = get(c,i) / (get(b,i) - get(a,i) * get(c,i-1));
			get(d,i) = (get(d,i) - get(d,i-1) * get(a,i)) / (get(b,i) - get(a,i) * get(c,i-1));  
		}

		get(x,num-1) = get(d,num-1);
	
		for (int i = num-2; i >= 0; i--)
			get(x,i) = get(d,i) - get(c,i) * get(x,i+1);
	}

	template<int dir, int var>
	__device__ void update_segment( FTYPE *x, Segment3D &seg, TimeLayer3D_GPU &layer, int id, int num_seg, int max_n )
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

			switch (dir)
			{
			case X: i++; break;
			case Y: j++; break;
			case Z: k++; break;
			case Z_as_Y: j++; break;
			}
		}
	}

	template<int dir, int var>
	__global__ void solve_segments( FluidParamsGPU p, int num_seg, Segment3D *segs, Node* nodes, TimeLayer3D_GPU cur, TimeLayer3D_GPU temp, TimeLayer3D_GPU next,
									FTYPE *a, FTYPE *b, FTYPE *c, FTYPE *d, FTYPE *x )
	{
		// fetch current segment info
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		if( id >= num_seg ) return;
		Segment3D &seg = segs[id];
		int max_n = max( max( cur.dimx, cur.dimy ), cur.dimz );
		int n = seg.size;

		apply_bc0<dir, var>(seg.posx, seg.posy, seg.posz, cur.dimy, cur.dimz, get(b,0), get(c,0), get(d,0), nodes);
		apply_bc1<dir, var>(seg.endx, seg.endy, seg.endz, cur.dimy, cur.dimz, get(a,n-1), get(b,n-1), get(d,n-1), nodes);
		
		build_matrix<dir, var>(p, seg.posx, seg.posy, seg.posz, a, b, c, d, n, cur, temp, id, num_seg, max_n);
			
		solve_tridiagonal(a, b, c, d, x, n, id, num_seg, max_n);
			
		update_segment<dir, var>(x, seg, next, id, num_seg, max_n);
	}

	template<int dir, int var>
	__global__ void build_matrix( FluidParamsGPU p, int num_seg, Segment3D *segs, Node* nodes, TimeLayer3D_GPU cur, TimeLayer3D_GPU temp, 
								  FTYPE *a, FTYPE *b, FTYPE *c, FTYPE *d )
	{
		// fetch current segment info
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		if( id >= num_seg ) return;
		Segment3D &seg = segs[id];
		int max_n = max( max( cur.dimx, cur.dimy ), cur.dimz );
		int n = seg.size;

		apply_bc0<dir, var>(seg.posx, seg.posy, seg.posz, cur.dimy, cur.dimz, get(b,0), get(c,0), get(d,0), nodes);
		apply_bc1<dir, var>(seg.endx, seg.endy, seg.endz, cur.dimy, cur.dimz, get(a,n-1), get(b,n-1), get(d,n-1), nodes);
		
		build_matrix<dir, var>(p, seg.posx, seg.posy, seg.posz, a, b, c, d, n, cur, temp, id, num_seg, max_n);	
	}

	template<int dir, int var>
	__global__ void solve_matrix( int num_seg, Segment3D *segs, TimeLayer3D_GPU next,
								  FTYPE *a, FTYPE *b, FTYPE *c, FTYPE *d, FTYPE *x )
	{
		// fetch current segment info
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		if( id >= num_seg ) return;
		Segment3D &seg = segs[id];
		int max_n = max( max( next.dimx, next.dimy ), next.dimz );
		int n = seg.size;

		solve_tridiagonal(a, b, c, d, x, n, id, num_seg, max_n);
			
		update_segment<dir, var>(x, seg, next, id, num_seg, max_n);
	}

	template<DirType dir, VarType var>
	void LaunchSolveSegments_dir_var( FluidParamsGPU p, int num_seg, Segment3D *segs, Node *nodes, TimeLayer3D_GPU &cur, TimeLayer3D_GPU &temp, TimeLayer3D_GPU &next,
									  FTYPE *d_a, FTYPE *d_b, FTYPE *d_c, FTYPE *d_d, FTYPE *d_x, bool decomposeOpt )
	{
		dim3 block(SOLVER_BLOCK_DIM);
		dim3 grid((num_seg + block.x - 1)/block.x);

		switch( decomposeOpt )
		{
		case true:
			build_matrix<dir, var><<<grid, block>>>( p, num_seg, segs, nodes, cur, temp, d_a, d_b, d_c, d_d );
			cudaThreadSynchronize();

			solve_matrix<dir, var><<<grid, block>>>( num_seg, segs, next, d_a, d_b, d_c, d_d, d_x );
			cudaThreadSynchronize();
			break;

		case false:
			cudaFuncSetCacheConfig(solve_segments<dir, var>, cudaFuncCachePreferL1);
			solve_segments<dir, var><<<grid, block>>>( p, num_seg, segs, nodes, cur, temp, next, d_a, d_b, d_c, d_d, d_x );
			cudaThreadSynchronize();
			break;
		}
	}

	template<DirType dir>
	void LaunchSolveSegments_dir( FluidParamsGPU p, int num_seg, Segment3D *segs, VarType var, Node *nodes, TimeLayer3D_GPU &cur, TimeLayer3D_GPU &temp, TimeLayer3D_GPU &next,
								  FTYPE *d_a, FTYPE *d_b, FTYPE *d_c, FTYPE *d_d, FTYPE *d_x, bool decomposeOpt )
	{
		switch( var )
		{
		case type_U: LaunchSolveSegments_dir_var<dir, type_U>( p, num_seg, segs, nodes, cur, temp, next, d_a, d_b, d_c, d_d, d_x, decomposeOpt ); break;
		case type_V: LaunchSolveSegments_dir_var<dir, type_V>( p, num_seg, segs, nodes, cur, temp, next, d_a, d_b, d_c, d_d, d_x, decomposeOpt ); break;
		case type_W: LaunchSolveSegments_dir_var<dir, type_W>( p, num_seg, segs, nodes, cur, temp, next, d_a, d_b, d_c, d_d, d_x, decomposeOpt ); break;
		case type_T: LaunchSolveSegments_dir_var<dir, type_T>( p, num_seg, segs, nodes, cur, temp, next, d_a, d_b, d_c, d_d, d_x, decomposeOpt ); break;
		}
	}

	void SolveSegments_GPU( FTYPE dt, FluidParams params, int num_seg, Segment3D *segs, VarType var, DirType dir, Node *nodes, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *next,
							FTYPE *d_a, FTYPE *d_b, FTYPE *d_c, FTYPE *d_d, FTYPE *d_x, bool decomposeOpt )
	{
		TimeLayer3D_GPU d_cur( cur );
		TimeLayer3D_GPU d_temp( temp );
		TimeLayer3D_GPU d_next( next );

		FluidParamsGPU p( var, dir, dt, cur->dx, cur->dy, cur->dz, params );

		switch( dir )
		{
		case X: LaunchSolveSegments_dir<X>( p, num_seg, segs, var, nodes, d_cur, d_temp, d_next, d_a, d_b, d_c, d_d, d_x, decomposeOpt ); break;
		case Y: LaunchSolveSegments_dir<Y>( p, num_seg, segs, var, nodes, d_cur, d_temp, d_next, d_a, d_b, d_c, d_d, d_x, decomposeOpt ); break;
		case Z: LaunchSolveSegments_dir<Z>( p, num_seg, segs, var, nodes, d_cur, d_temp, d_next, d_a, d_b, d_c, d_d, d_x, decomposeOpt ); break;
		case Z_as_Y: LaunchSolveSegments_dir<Z_as_Y>( p, num_seg, segs, var, nodes, d_cur, d_temp, d_next, d_a, d_b, d_c, d_d, d_x, decomposeOpt ); break;
		}
	}
}