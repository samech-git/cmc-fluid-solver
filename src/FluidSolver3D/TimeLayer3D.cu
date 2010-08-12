// GPU implementation of TimeLayer3D

#include "Grid3D.h"

#define BLOCK_DIM_X		32
#define BLOCK_DIM_Y		8

using namespace FluidSolver3D;

__global__ void copy_grid(int dimx, int dimy, int dimz, FTYPE *u, FTYPE *v, FTYPE *w, FTYPE *T, Node *nodes, NodeType target)
{
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;

	if( k >= dimz || j >= dimy ) return;

	for (int i = 0; i < dimx; i++)
	{
		int id = i * dimy * dimz + j * dimz + k;
		if( nodes[id].type == target )
		{
			u[id] = nodes[id].v.x;
			v[id] = nodes[id].v.y;
			w[id] = nodes[id].v.z;
			T[id] = nodes[id].T;
		}
	}
}

__global__ void clear(int dimx, int dimy, int dimz, FTYPE *u, FTYPE *v, FTYPE *w, FTYPE *T, Node *nodes, NodeType target, FTYPE const_u, FTYPE const_v, FTYPE const_w, FTYPE const_T)
{
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;

	if( k >= dimz || j >= dimy ) return;

	for (int i = 0; i < dimx; i++)
	{
		int id = i * dimy * dimz + j * dimz + k;
		if( nodes[id].type == target )
		{
			u[id] = const_u;
			v[id] = const_v;
			w[id] = const_w;
			T[id] = const_T;
		}
	}
}

__global__ void copy(int dimx, int dimy, int dimz, FTYPE *src, FTYPE *dest, Node *nodes, NodeType target)
{
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;

	if( k >= dimz || j >= dimy ) return;

	for (int i = 0; i < dimx; i++)
	{
		int id = i * dimy * dimz + j * dimz + k;
		if( nodes[id].type == target )
			dest[id] = src[id];
	}
}

__global__ void merge(int dimx, int dimy, int dimz, FTYPE *src, FTYPE *dest, Node *nodes, NodeType target)
{
	int k = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;

	if( k >= dimz || j >= dimy ) return;

	for (int i = 0; i < dimx; i++)
	{
		int id = i * dimy * dimz + j * dimz + k;
		if( nodes[id].type == target )
			dest[id] = ( dest[id] + src[id] ) / 2;
	}
}

void CopyFromGrid_GPU(int dimx, int dimy, int dimz, FTYPE *u, FTYPE *v, FTYPE *w, FTYPE *T, Node *nodes, NodeType target)
{
	dim3 block(BLOCK_DIM_X, BLOCK_DIM_Y);
	dim3 grid((dimz + block.x - 1)/block.x, (dimy + block.y - 1)/block.y);
	copy_grid<<<grid, block>>>(dimx, dimy, dimz, u, v, w, T, nodes, target);
	cudaThreadSynchronize();
}

void CopyFieldTo_GPU(int dimx, int dimy, int dimz, FTYPE *src, FTYPE *dest, Node *nodes, NodeType target)
{
	dim3 block(BLOCK_DIM_X, BLOCK_DIM_Y);
	dim3 grid((dimz + block.x - 1)/block.x, (dimy + block.y - 1)/block.y);
	copy<<<grid, block>>>(dimx, dimy, dimz, src, dest, nodes, target);
	cudaThreadSynchronize();
}

void MergeFieldTo_GPU(int dimx, int dimy, int dimz, FTYPE *src, FTYPE *dest, Node *nodes, NodeType target)
{
	dim3 block(BLOCK_DIM_X, BLOCK_DIM_Y);
	dim3 grid((dimz + block.x - 1)/block.x, (dimy + block.y - 1)/block.y);
	merge<<<grid, block>>>(dimx, dimy, dimz, src, dest, nodes, target);
	cudaThreadSynchronize();
}

void Clear_GPU(int dimx, int dimy, int dimz, FTYPE *u, FTYPE *v, FTYPE *w, FTYPE *T, Node *nodes, NodeType target, FTYPE const_u, FTYPE const_v, FTYPE const_w, FTYPE const_T)
{
	dim3 block(BLOCK_DIM_X, BLOCK_DIM_Y);
	dim3 grid((dimz + block.x - 1)/block.x, (dimy + block.y - 1)/block.y);
	clear<<<grid, block>>>(dimx, dimy, dimz, u, v, w, T, nodes, target, const_u, const_v, const_w, const_T);
	cudaThreadSynchronize();
}