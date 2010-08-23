// GPU implementation of TimeLayer3D

#include "Grid3D.h"

#define COPY_BLOCK_DIM_X		32
#define COPY_BLOCK_DIM_Y		8

#define TRANSPOSE_TILE_DIM		16
#define TRANSPOSE_BLOCK_ROWS	16

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

__global__ void transpose(int dimx, int dimy, int dimz, FTYPE *src, FTYPE *dest)
{
	__shared__ FTYPE tile[TRANSPOSE_TILE_DIM][TRANSPOSE_TILE_DIM+1];

    // read the tile from global memory into shared memory
	int k0 = blockIdx.x * TRANSPOSE_TILE_DIM + threadIdx.x;
	int j0 = blockIdx.y * TRANSPOSE_TILE_DIM + threadIdx.y;

	int j1 = blockIdx.y * TRANSPOSE_TILE_DIM + threadIdx.x;
    int k1 = blockIdx.x * TRANSPOSE_TILE_DIM + threadIdx.y;

	int base_u0 = k0 + j0 * dimz;	
    int base_u1 = j1 + k1 * dimy;
	
	for (int i = 0; i < dimx; i++)
	{	
		// read tile from global to shared memory
		for (int row = 0; row < TRANSPOSE_TILE_DIM; row += TRANSPOSE_BLOCK_ROWS)
			tile[threadIdx.y + row][threadIdx.x] = src[base_u0 + row * dimz];
		base_u0 += dimz * dimy;
	
		__syncthreads();
	
		// write the transposed tile to global memory 
		for (int row = 0; row < TRANSPOSE_TILE_DIM; row += TRANSPOSE_BLOCK_ROWS)
			dest[base_u1 + row * dimy] = tile[threadIdx.x][threadIdx.y + row];
		base_u1 += dimz * dimy;
	}
}


void CopyFromGrid_GPU(int dimx, int dimy, int dimz, FTYPE *u, FTYPE *v, FTYPE *w, FTYPE *T, Node *nodes, NodeType target)
{
	dim3 block(COPY_BLOCK_DIM_X, COPY_BLOCK_DIM_Y);
	dim3 grid((dimz + block.x - 1)/block.x, (dimy + block.y - 1)/block.y);
	copy_grid<<<grid, block>>>(dimx, dimy, dimz, u, v, w, T, nodes, target);
	cudaThreadSynchronize();
}

void CopyFieldTo_GPU(int dimx, int dimy, int dimz, FTYPE *src, FTYPE *dest, Node *nodes, NodeType target)
{
	dim3 block(COPY_BLOCK_DIM_X, COPY_BLOCK_DIM_Y);
	dim3 grid((dimz + block.x - 1)/block.x, (dimy + block.y - 1)/block.y);
	copy<<<grid, block>>>(dimx, dimy, dimz, src, dest, nodes, target);
	cudaThreadSynchronize();
}

void MergeFieldTo_GPU(int dimx, int dimy, int dimz, FTYPE *src, FTYPE *dest, Node *nodes, NodeType target)
{
	dim3 block(COPY_BLOCK_DIM_X, COPY_BLOCK_DIM_Y);
	dim3 grid((dimz + block.x - 1)/block.x, (dimy + block.y - 1)/block.y);
	merge<<<grid, block>>>(dimx, dimy, dimz, src, dest, nodes, target);
	cudaThreadSynchronize();
}

void Clear_GPU(int dimx, int dimy, int dimz, FTYPE *u, FTYPE *v, FTYPE *w, FTYPE *T, Node *nodes, NodeType target, FTYPE const_u, FTYPE const_v, FTYPE const_w, FTYPE const_T)
{
	dim3 block(COPY_BLOCK_DIM_X, COPY_BLOCK_DIM_Y);
	dim3 grid((dimz + block.x - 1)/block.x, (dimy + block.y - 1)/block.y);
	clear<<<grid, block>>>(dimx, dimy, dimz, u, v, w, T, nodes, target, const_u, const_v, const_w, const_T);
	cudaThreadSynchronize();
}

void Transpose_GPU(int dimx, int dimy, int dimz, FTYPE *u, FTYPE *dest_u)
{
	dim3 block(TRANSPOSE_TILE_DIM, TRANSPOSE_BLOCK_ROWS);
	dim3 grid((dimz + TRANSPOSE_TILE_DIM - 1)/TRANSPOSE_TILE_DIM, (dimy + TRANSPOSE_TILE_DIM - 1)/TRANSPOSE_TILE_DIM);
	transpose<<<grid, block>>>(dimx, dimy, dimz, u, dest_u);
	cudaThreadSynchronize();
}