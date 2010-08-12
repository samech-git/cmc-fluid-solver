#include "AdiSolver3D.h"

#define BLOCK_DIM		128

namespace FluidSolver3D
{
	struct FluidParamsGPU
	{
		FTYPE vis_dx2, vis_dy2, vis_dz2;
		FTYPE dt, dx, dy, dz;
		FTYPE v_T, t_phi;

		FluidParamsGPU( VarType var, FTYPE _dt, FTYPE _dx, FTYPE _dy, FTYPE _dz, FluidParams _params ) :
			dt(_dt), dx(_dx), dy(_dy), dz(_dz)
		{
			switch (var)
			{
			case type_U:
			case type_V:
			case type_W:
				vis_dx2 = _params.v_vis / (dx * dx);
				vis_dy2 = _params.v_vis / (dy * dy);
				vis_dz2 = _params.v_vis / (dz * dz);
				break;
			case type_T:
				vis_dx2 = _params.t_vis / (dx * dx);
				vis_dy2 = _params.t_vis / (dy * dy);
				vis_dz2 = _params.t_vis / (dz * dz);
				break;
			}
						
			v_T = _params.v_T;
			t_phi = _params.t_phi;
		}
	};

	/*#define elem(u, i, j, k)	(u[i * dimy * dimz + j * dimz + k])
	
	#define d_x(u, i, j, k)		((elem(u, i+1, j, k) - elem(u, i-1, j, k)) / (2 * dx))
	#define d_y(u, i, j, k)		((elem(u, i, j+1, k) - elem(u, i, j-1, k)) / (2 * dy))
	#define d_z(u, i, j, k)		((elem(u, i, j, k+1) - elem(u, i, j, k-1)) / (2 * dz))
	#define d_xx(u, i, j, k)	((elem(u, i+1, j, k) - 2 * elem(u, i, j, k) + elem(u, i-1, j, k)) / (dx * dx))
	#define d_yy(u, i, j, k)	((elem(u, i, j+1, k) - 2 * elem(u, i, j, k) + elem(u, i, j-1, k)) / (dy * dy))
	#define d_zz(u, i, j, k)	((elem(u, i, j, k+1) - 2 * elem(u, i, j, k) + elem(u, i, j, k-1)) / (dz * dz))*/

	//#define DissFuncX(u, v, w, i, j, k) (2 * d_x(u, i, j, k) * d_x(u, i, j, k) + d_x(v, i, j, k) * d_x(v, i, j, k) + d_x(w, i, j, k) * d_x(w, i, j, k) + d_x(v, i, j, k) * d_y(u, i, j, k) + d_x(w, i, j, k) * d_z(u, i, j, k))

	template<int var>
	__device__ void apply_bc0(int i, int j, int k, int dimy, int dimz, FTYPE &b0, FTYPE &c0, FTYPE &d0, Node *nodes)
	{
		int id = i * dimy * dimz + j * dimz + k;
		if ((var == type_T && nodes[id].bc_temp == BC_FREE) ||
			(var != type_T && nodes[id].bc_vel == BC_FREE))
		{
			// free
			b0 = 1.0; 
			c0 = -1.0; 
			d0 = 0.0; 
		}
		else
		{
			// no-slip
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

	template<int var>
	__device__ void apply_bc1(int i, int j, int k, int dimy, int dimz, FTYPE &a1, FTYPE &b1, FTYPE &d1, Node *nodes)
	{
		int id = i * dimy * dimz + j * dimz + k;
		if ((var == type_T && nodes[id].bc_temp == BC_FREE) ||
			(var != type_T && nodes[id].bc_vel == BC_FREE))
		{
			// free
			a1 = 1.0; 
			b1 = -1.0; 
			d1 = 0.0;
		}
		else
		{
			// no-slip
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
	__device__ void build_matrix(FluidParamsGPU params, int i, int j, int k, FTYPE *a, FTYPE *b, FTYPE *c, FTYPE *d, int n, TimeLayer3D_GPU &cur, TimeLayer3D_GPU &temp)
	{	
		for (int p = 1; p < n-1; p++)
		{
			switch (dir)
			{
			case X:		
				a[p] = - temp.elem(temp.u, i+p, j, k) / (2 * temp.dx) - params.vis_dx2; 
				b[p] = 3 / params.dt  +  2 * params.vis_dx2; 
				c[p] = temp.elem(temp.u, i+p, j, k) / (2 * temp.dx) - params.vis_dx2; 
				
				switch (var)	
				{
				case type_U: d[p] = cur.elem(cur.u, i+p, j, k) * 3 / params.dt - params.v_T * temp.d_x(temp.T, i+p, j, k); break;
				case type_V: d[p] = cur.elem(cur.v, i+p, j, k) * 3 / params.dt; break;
				case type_W: d[p] = cur.elem(cur.w, i+p, j, k) * 3 / params.dt; break;
				case type_T: d[p] = cur.elem(cur.T, i+p, j, k) * 3 / params.dt + params.t_phi * temp.DissFuncX(i+p, j, k); break;
				}	
				break;

			case Y:
				a[p] = - temp.elem(temp.v, i, j+p, k) / (2 * temp.dy) - params.vis_dy2; 
				b[p] = 3 / params.dt  +  2 * params.vis_dy2; 
				c[p] = temp.elem(temp.v, i, j+p, k) / (2 * temp.dy) - params.vis_dy2; 
				
				switch (var)	
				{
				case type_U: d[p] = cur.elem(cur.u, i, j+p, k) * 3 / params.dt; break;
				case type_V: d[p] = cur.elem(cur.v, i, j+p, k) * 3 / params.dt - params.v_T * temp.d_y(temp.T, i, j+p, k); break;
				case type_W: d[p] = cur.elem(cur.w, i, j+p, k) * 3 / params.dt; break;
				case type_T: d[p] = cur.elem(cur.T, i, j+p, k) * 3 / params.dt + params.t_phi * temp.DissFuncY(i, j+p, k); break;
				}
				break;

			case Z:
				a[p] = - temp.elem(temp.w, i, j, k+p) / (2 * temp.dz) - params.vis_dz2; 
				b[p] = 3 / params.dt  +  2 * params.vis_dz2; 
				c[p] = temp.elem(temp.w, i, j, k+p) / (2 * temp.dz) - params.vis_dz2; 
				
				switch (var)	
				{
				case type_U: d[p] = cur.elem(cur.u, i, j, k+p) * 3 / params.dt; break;
				case type_V: d[p] = cur.elem(cur.v, i, j, k+p) * 3 / params.dt; break;
				case type_W: d[p] = cur.elem(cur.w, i, j, k+p) * 3 / params.dt - params.v_T * temp.d_z(temp.T, i, j, k+p); break;
				case type_T: d[p] = cur.elem(cur.T, i, j, k+p) * 3 / params.dt + params.t_phi * temp.DissFuncZ(i, j, k+p); break;
				}
				break;
			}
		}
	}

	__device__ void solve_tridiagonal( FTYPE *a, FTYPE *b, FTYPE *c, FTYPE *d, FTYPE *x, int num )
	{
		c[num-1] = 0.0;
		
		c[0] = c[0] / b[0];
		d[0] = d[0] / b[0];

		for (int i = 1; i < num; i++)
		{
			c[i] = c[i] / (b[i] - a[i] * c[i-1]);
			d[i] = (d[i] - d[i-1] * a[i]) / (b[i] - a[i] * c[i-1]);  
		}

		x[num-1] = d[num-1];
	
		for (int i = num-2; i >= 0; i--)
			x[i] = d[i] - c[i] * x[i+1];
	}

	template<int dir, int var>
	__device__ void update_segment( FTYPE *d_x, Segment3D &seg, TimeLayer3D_GPU &layer )
	{
		int i = seg.posx;
		int j = seg.posy;
		int k = seg.posz;

		for (int t = 0; t < seg.size; t++)
		{
			switch (var)
			{
			case type_U: layer.elem(layer.u, i, j, k) = d_x[t]; break;
			case type_V: layer.elem(layer.v, i, j, k) = d_x[t]; break;
			case type_W: layer.elem(layer.w, i, j, k) = d_x[t]; break;
			case type_T: layer.elem(layer.T, i, j, k) = d_x[t]; break;
			}

			switch (dir)
			{
			case X: i++; break;
			case Y: j++; break;
			case Z: k++; break;
			}
		}
	}

	template<int dir, int var>
	__global__ void solve_segments( FluidParamsGPU p, int num_seg, Segment3D *segs, Node* nodes, TimeLayer3D_GPU cur, TimeLayer3D_GPU temp, TimeLayer3D_GPU next,
									FTYPE *d_a, FTYPE *d_b, FTYPE *d_c, FTYPE *d_d, FTYPE *d_x )
	{
		// fetch current segment info
		int id = blockIdx.x * blockDim.x + threadIdx.x;
		if( id > num_seg ) return;
		Segment3D &seg = segs[id];
		int max_n = max(cur.dimx, max(cur.dimy, cur.dimz));
		int n = seg.size;

		FTYPE *a = d_a + id * max_n;
		FTYPE *b = d_b + id * max_n;
		FTYPE *c = d_c + id * max_n;
		FTYPE *d = d_d + id * max_n;
		FTYPE *x = d_x + id * max_n;

		apply_bc0<var>(seg.posx, seg.posy, seg.posz, cur.dimy, cur.dimz, b[0], c[0], d[0], nodes);
		apply_bc1<var>(seg.endx, seg.endy, seg.endz, cur.dimy, cur.dimz, a[n-1], b[n-1], d[n-1], nodes);
		
		build_matrix<dir, var>(p, seg.posx, seg.posy, seg.posz, a, b, c, d, n, cur, temp);
			
		solve_tridiagonal(a, b, c, d, x, n);
			
		update_segment<dir, var>(x, seg, next);
	}

	template<DirType dir, VarType var>
	void LaunchSolveSegments_dir_var( FluidParamsGPU p, int num_seg, Segment3D *segs, Node* nodes, TimeLayer3D_GPU &cur, TimeLayer3D_GPU &temp, TimeLayer3D_GPU &next,
									  FTYPE *d_a, FTYPE *d_b, FTYPE *d_c, FTYPE *d_d, FTYPE *d_x )
	{
		dim3 block(BLOCK_DIM);
		dim3 grid((num_seg + block.x - 1)/block.x);
		solve_segments<dir, var><<<grid, block>>>( p, num_seg, segs, nodes, cur, temp, next, d_a, d_b, d_c, d_d, d_x );
		cudaThreadSynchronize();
	}

	template<DirType dir>
	void LaunchSolveSegments_dir( FluidParamsGPU p, int num_seg, Segment3D *segs, VarType var, Node* nodes, TimeLayer3D_GPU &cur, TimeLayer3D_GPU &temp, TimeLayer3D_GPU &next,
								  FTYPE *d_a, FTYPE *d_b, FTYPE *d_c, FTYPE *d_d, FTYPE *d_x )
	{
		switch( var )
		{
		case type_U: LaunchSolveSegments_dir_var<dir, type_U>( p, num_seg, segs, nodes, cur, temp, next, d_a, d_b, d_c, d_d, d_x ); break;
		case type_V: LaunchSolveSegments_dir_var<dir, type_V>( p, num_seg, segs, nodes, cur, temp, next, d_a, d_b, d_c, d_d, d_x ); break;
		case type_W: LaunchSolveSegments_dir_var<dir, type_W>( p, num_seg, segs, nodes, cur, temp, next, d_a, d_b, d_c, d_d, d_x ); break;
		case type_T: LaunchSolveSegments_dir_var<dir, type_T>( p, num_seg, segs, nodes, cur, temp, next, d_a, d_b, d_c, d_d, d_x ); break;
		}
	}

	void SolveSegments_GPU( FTYPE dt, FluidParams params, int num_seg, Segment3D *segs, VarType var, DirType dir, Node* nodes, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *next,
							FTYPE *d_a, FTYPE *d_b, FTYPE *d_c, FTYPE *d_d, FTYPE *d_x )
	{
		TimeLayer3D_GPU d_cur( cur );
		TimeLayer3D_GPU d_temp( temp );
		TimeLayer3D_GPU d_next( next );

		FluidParamsGPU p( var, dt, cur->dx, cur->dy, cur->dz, params );

		switch( dir )
		{
		case X: LaunchSolveSegments_dir<X>( p, num_seg, segs, var, nodes, d_cur, d_temp, d_next, d_a, d_b, d_c, d_d, d_x ); break;
		case Y: LaunchSolveSegments_dir<Y>( p, num_seg, segs, var, nodes, d_cur, d_temp, d_next, d_a, d_b, d_c, d_d, d_x ); break;
		case Z: LaunchSolveSegments_dir<Z>( p, num_seg, segs, var, nodes, d_cur, d_temp, d_next, d_a, d_b, d_c, d_d, d_x ); break;
		}
	}
}