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

#ifdef linux
#include <cmath>  // for abs functions
#endif

using namespace FluidSolver3D;

#include <cuda_runtime.h>

extern void CopyFieldTo_GPU(int dimx, int dimy, int dimz, FTYPE **src, FTYPE **dest, NodeType **nodes, NodeType target, int haloSize);
extern void MergeFieldTo_GPU(int dimx, int dimy, int dimz, FTYPE **src, FTYPE **dest, NodeType **nodes, NodeType target, int haloSize);
extern void CopyFromGrid_GPU(int dimx, int dimy, int dimz, FTYPE **u, FTYPE **v, FTYPE **w, FTYPE **T, Node **nodes, NodeType target, int haloSize);
extern void CopyGridBoundary_GPU(DirType dir, int dimx, int dimy, int dimz, FTYPE **u, FTYPE **v, FTYPE **w, FTYPE **T, int *num_seg, Segment3D **segs, NodesBoundary3D **nodes, int haloSize);
extern void Clear_GPU(int dimx, int dimy, int dimz, FTYPE **u, FTYPE **v, FTYPE **w, FTYPE **T, NodeType **nodes, NodeType target, FTYPE const_u, FTYPE const_v, FTYPE const_w, FTYPE const_T, int haloSize);
extern void Transpose_GPU(int dimx, int dimy, int dimz, FTYPE **u, FTYPE **dest_u, int haloSize); // need to test

namespace FluidSolver3D
{

	enum SwipeType
	{
		ALL,
		FORWARD,
		BACK
	};

#ifdef __PARA
template <typename T, SwipeType swipe>
	void paraSend(T* src, int num_elems, int tagID = 666, MPI_Request *request = NULL)
	{
		PARAplan* pplan = PARAplan::Instance();
		int irank =  pplan->rank();
		int size = pplan->size();
		switch (swipe)
		{
		case FORWARD:
			if (irank < size - 1)
				if (request == NULL)
					mpiSafeCall(MPI_Send(src, num_elems, mpi_typeof(src), irank + 1, tagID, MPI_COMM_WORLD), "paraSend<FORWARD>: MPI_Send");
				else
					mpiSafeCall(MPI_Isend(src, num_elems, mpi_typeof(src), irank + 1, tagID, MPI_COMM_WORLD, request), "paraSend<FORWARD>: MPI_Isend");
			break;
		case BACK:
			if (irank > 0)
				if (request == NULL)
					mpiSafeCall(MPI_Send(src, num_elems, mpi_typeof(src), irank - 1, tagID, MPI_COMM_WORLD), "paraSend<BACK>: MPI_Send");
				else
					mpiSafeCall(MPI_Isend(src, num_elems, mpi_typeof(src), irank - 1, tagID, MPI_COMM_WORLD, request), "paraSend<BACK>: MPI_Isend");
			break;
		}
	}
#endif
	
#ifdef __PARA
template <typename T, SwipeType swipe>
	void paraRecv(T* dst, int num_elems, int tagID = 666, MPI_Request *request = NULL)
	{
		PARAplan* pplan = PARAplan::Instance();
		int irank =  pplan->rank();
		int size = pplan->size();
		MPI_Status status;
		switch (swipe)
		{
		case FORWARD:
			if (irank > 0)
				if (request == NULL)
					mpiSafeCall(MPI_Recv(dst, num_elems, mpi_typeof(dst), irank-1, tagID, MPI_COMM_WORLD, &status), "paraRecv<FORWARD>: MPI_Recv");
				else
					mpiSafeCall(MPI_Irecv(dst, num_elems, mpi_typeof(dst), irank-1, tagID, MPI_COMM_WORLD, request), "paraRecv<FORWARD>: MPI_Irecv");
			break;
		case BACK:
			if (irank < size - 1)
				if (request == NULL)
					mpiSafeCall(MPI_Recv(dst, num_elems, mpi_typeof(dst), irank+1, tagID, MPI_COMM_WORLD, &status), "paraRecv<BACK>: MPI_Recv");
				else
					mpiSafeCall(MPI_Irecv(dst, num_elems, mpi_typeof(dst), irank+1, tagID, MPI_COMM_WORLD, request), "paraRecv<BACK>: MPI_Irecv");
			break;
		}
	}
#endif

template <typename T, SwipeType swipe>
	void paraDevSend(T *dev_src, T *mpi_buf, int num_elems, int tagID = 666)
	{		
#ifdef __PARA
		PARAplan* pplan = PARAplan::Instance();
		int irank =  pplan->rank();
		int size = pplan->size();
		switch (swipe)
		{
		case FORWARD:
			if (irank < size - 1)
			{
				gpuSafeCall(cudaMemcpy(mpi_buf, dev_src, sizeof(T) * num_elems, cudaMemcpyDeviceToHost), "paraDevSend: cudaMemcpy");
				cudaDeviceSynchronize();
				mpiSafeCall(MPI_Ssend(mpi_buf, num_elems, mpi_typeof(mpi_buf), irank + 1, tagID, MPI_COMM_WORLD), "paraDevSend: MPI_Ssend");
			}
			break;
		case BACK:
			if (irank > 0)
			{
				cudaMemcpy(mpi_buf, dev_src, sizeof(T) * num_elems, cudaMemcpyDeviceToHost);
				cudaDeviceSynchronize();
				mpiSafeCall(MPI_Ssend(mpi_buf, num_elems, mpi_typeof(mpi_buf), irank - 1, tagID, MPI_COMM_WORLD), "paraDevSend: MPI_Ssend");
			}
			break;
		}
#endif
	}

template <typename T, SwipeType swipe>
	void paraDevRecv(T* dev_dst, T *mpi_buf, int num_elems, int tagID = 666)
	{
#ifdef __PARA
		PARAplan* pplan = PARAplan::Instance();
		int irank =  pplan->rank();
		int size = pplan->size();
		MPI_Status status;
		switch (swipe)
		{
		case FORWARD:
			if (irank > 0)
			{
				mpiSafeCall(MPI_Recv(mpi_buf, num_elems, mpi_typeof(mpi_buf), irank-1, tagID, MPI_COMM_WORLD, &status),  "paraDevRecv: MPI_Recv");
				cudaMemcpy(dev_dst, mpi_buf, sizeof(T) * num_elems, cudaMemcpyHostToDevice);
			}
			break;
		case BACK:
			if (irank < size - 1)
			{
				MPI_Recv(mpi_buf, num_elems, mpi_typeof(mpi_buf), irank+1, tagID, MPI_COMM_WORLD, &status);
				cudaMemcpy(dev_dst, mpi_buf, sizeof(T) * num_elems, cudaMemcpyHostToDevice);
			}
			break;
		}
#endif
	}

template <typename T, SwipeType swipe>
	void haloMemcpyPeer(T** dev_array, int idev, int haloSize, int num_elems_total, int haloShift = 0, int haloPartSize = -1)
	/*
		propagates halo of multi device array back or forth
		num_elems_total - number of elements without halo
		haloPartSize and haloShift are required to sent fractions of the halo
	*/
	{
		int ip1 = idev+1; int i = idev; int im1 = idev-1;
		if (haloPartSize == -1)
			haloPartSize = haloSize;
	#if MGPU_EMU
		im1 = ip1 = i = DEFAULT_DEVICE;
	#endif
		switch (swipe)
		{
		case FORWARD:
			gpuSafeCall(cudaMemcpyPeer(dev_array[idev+1] + haloShift,  ip1, dev_array[idev] + haloSize + (num_elems_total - haloSize) + haloShift, i,  sizeof(T) * haloPartSize), "haloMemcpyPeer: FORWARD", i);
			break;
		case BACK:
			gpuSafeCall(cudaMemcpyPeer(dev_array[idev-1] + haloSize + num_elems_total + haloShift,  im1, dev_array[idev] + haloSize + haloShift, i,  sizeof(T) * haloPartSize), "haloMemcpyPeer: BACK", i);
			break;
		}
	}

template <typename T, SwipeType swipe>
	void haloMemcpyPeerAsync(T** dev_array, int idev, int haloSize, int num_elems_total, cudaStream_t& stream, int haloShift = 0, int haloPartSize = -1)
	/*
		propagates halo of multi device array back or forth
		num_elems_total - number of elements without halo
	*/
	{
		int ip1 = idev+1; int i = idev; int im1 = idev-1;
		if (haloPartSize == -1)
			haloPartSize = haloSize;
	#if MGPU_EMU
		im1 = ip1 = i = DEFAULT_DEVICE;
	#endif
		switch (swipe)
		{
		case FORWARD:
			gpuSafeCall( cudaMemcpyPeerAsync(dev_array[idev+1] + haloShift,  ip1, dev_array[idev] + haloSize + (num_elems_total - haloSize) + haloShift, i,  sizeof(T) * haloPartSize, stream), "haloMemcpyPeerAsync: FORWARD", i);
			break;
		case BACK:
			gpuSafeCall( cudaMemcpyPeerAsync(dev_array[idev-1] + haloSize + num_elems_total + haloShift,  im1, dev_array[idev] + haloSize + haloShift, i,  sizeof(T) * haloPartSize, stream), "haloMemcpyPeerAsync: BACK", i);
			break;
		}
	}  

	struct ScalarField3D
	{
		BackendType hw; 
		int dimx, dimy, dimz;
		FTYPE dx, dy, dz;

		// access element
		inline FTYPE& elem(int i, int j, int k)
		{
			return u[haloSize + i * dimy * dimz + j * dimz + k];
		}

		// access the whole array
		FTYPE *getArray()
		{
			return u;
		}

		FTYPE **getMultiArray()
		{
			return  dd_u;
		}

		void syncHalos(int tagID_F = 666, int tagID_B = 667, FTYPE* mpi_buf = NULL)
		/*
			Will synchronize fields' halos if haloSize != 0 
		*/
		{
			if (haloSize <= 0)
				return;
			PARAplan *pplan = PARAplan::Instance();
			switch (hw)
			{
			case CPU:
#ifdef __PARA
				if (pplan->size() > 1)
				{
					paraRecv<FTYPE, FORWARD>(u, haloSize, tagID_F);
					paraSend<FTYPE, FORWARD>(u + haloSize + pplan->getLength1D() * haloSize - haloSize, haloSize, tagID_F);

					paraRecv<FTYPE, BACK>(u + haloSize +  pplan->getLength1D() * haloSize, haloSize, tagID_B);
					paraSend<FTYPE, BACK>(u + haloSize, haloSize, tagID_B);
				}
#endif
				break;
			case GPU:
				GPUplan *pGPUplan = GPUplan::Instance();
				int gpuSize = pGPUplan->size();
				for (int i = 0; i < gpuSize; i++)
				{
					if ( i < gpuSize - 1)
						haloMemcpyPeerAsync<FTYPE, FORWARD>( dd_u, i, haloSize, haloSize * pGPUplan->node(i)->getLength1D(), pGPUplan->node(i)->stream);
					if ( i > 0)
						haloMemcpyPeerAsync<FTYPE, BACK>( dd_u, i, haloSize, haloSize * pGPUplan->node(i-1)->getLength1D(), pGPUplan->node(i)->stream2);
				}
				for (int i = 0; i < gpuSize; i++)
				{
					pGPUplan->setDevice(i);
					gpuSafeCall(cudaStreamSynchronize(pGPUplan->node(i)->stream), "ScalarField3D::syncHalos(): cudaStreamSynchronize", i);
					gpuSafeCall(cudaStreamSynchronize(pGPUplan->node(i)->stream2), "ScalarField3D::syncHalos(): cudaStreamSynchronize", i);
				}
#ifdef __PARA				
				int irank = pplan->rank();
				int size = pplan->size();
				if (size > 1)
				{
					FTYPE *mpi_buf_s = mpi_buf;
					FTYPE *mpi_buf_r = mpi_buf + haloSize;
					MPI_Request request_s, request_r;
					MPI_Status status;
					pGPUplan->setDevice(gpuSize-1);
					gpuSafeCall( cudaMemcpy(mpi_buf_s, dd_u[gpuSize-1] + haloSize + pGPUplan->node(gpuSize - 1)->getLength1D() * haloSize - haloSize, sizeof(FTYPE) * haloSize, cudaMemcpyDeviceToHost), "syncHalos<FORWARD>: cudaMemcpy" );
					if (irank > 0)
						mpiSafeCall( MPI_Irecv(mpi_buf_r, haloSize, mpi_typeof(mpi_buf_r), irank - 1, tagID_F, MPI_COMM_WORLD, &request_r), "syncHalos<FORWARD>: MPI_Irecv" );
					if (irank < size - 1)
						mpiSafeCall( MPI_Isend(mpi_buf_s, haloSize, mpi_typeof(mpi_buf_s), irank + 1, tagID_F, MPI_COMM_WORLD, &request_s), "syncHalos<FORWARD>: MPI_Isend" );
					if (irank > 0)
						MPI_Wait(&request_r, &status);
					if (irank < size - 1)
						MPI_Wait(&request_s, &status);				
					pGPUplan->setDevice(0);
					gpuSafeCall( cudaMemcpy(dd_u[0], mpi_buf_r, sizeof(FTYPE) * haloSize, cudaMemcpyHostToDevice), "syncHalos<FORWARD>: cudaMemcpy" );
					gpuSafeCall( cudaMemcpy(mpi_buf_s, dd_u[0] + haloSize, sizeof(FTYPE) * haloSize, cudaMemcpyDeviceToHost), "syncHalos<BACK>: cudaMemcpy" );
					if (irank < size - 1)
						mpiSafeCall( MPI_Irecv(mpi_buf_r, haloSize, mpi_typeof(mpi_buf_r), irank + 1, tagID_B, MPI_COMM_WORLD, &request_r), "syncHalos<BACK>: MPI_Irecv" );
					if (irank > 0)
						mpiSafeCall( MPI_Isend(mpi_buf_s, haloSize, mpi_typeof(mpi_buf_s), irank - 1, tagID_B, MPI_COMM_WORLD, &request_s), "syncHalos<BACK>: MPI_Isend" );
					if (irank < size - 1)
						MPI_Wait(&request_r, &status);
					if (irank > 0)
						MPI_Wait(&request_s, &status);	
					pGPUplan->setDevice(gpuSize-1);
					gpuSafeCall( cudaMemcpy(dd_u[gpuSize-1] + haloSize +  pGPUplan->node(gpuSize - 1)->getLength1D() * haloSize, mpi_buf_r, sizeof(FTYPE) * haloSize, cudaMemcpyHostToDevice), "syncHalos<BACK>: cudaMemcpy" );
				}
#endif
				pGPUplan->deviceSynchronize();
				break;
			}
		}

		// get derivatives
		inline FTYPE d_x(int i, int j, int k)	{ return (elem(i+1, j, k) - elem(i-1, j, k)) / (2 * dx); }
		inline FTYPE d_y(int i, int j, int k)	{ return (elem(i, j+1, k) - elem(i, j-1, k)) / (2 * dy); }
		inline FTYPE d_z(int i, int j, int k)	{ return (elem(i, j, k+1) - elem(i, j, k-1)) / (2 * dz); }
		inline FTYPE d_xx(int i, int j, int k)	{ return (elem(i+1, j, k) - 2 * elem(i, j, k) + elem(i-1, j, k)) / (dx * dx); }
		inline FTYPE d_yy(int i, int j, int k)	{ return (elem(i, j+1, k) - 2 * elem(i, j, k) + elem(i, j-1, k)) / (dy * dy); }
		inline FTYPE d_zz(int i, int j, int k)	{ return (elem(i, j, k+1) - 2 * elem(i, j, k) + elem(i, j, k-1)) / (dz * dz); }

		ScalarField3D(BackendType _hw, int _dimx, int _dimy, int _dimz, FTYPE _dx, FTYPE _dy, FTYPE _dz, int _haloSize = 0) : 
			hw(_hw), dimx(_dimx), dimy(_dimy), dimz(_dimz),
			dx(_dx), dy(_dy), dz(_dz), haloSize(_haloSize)
		{
			PARAplan *pplan = PARAplan::Instance();
			dimxOffset = pplan->getOffset1D();
			switch( hw )
			{
			case CPU: u = new FTYPE[dimx * dimy * dimz + 2 * haloSize]; break;
			case GPU: multiDevAlloc<FTYPE>(dd_u, dimx * dimy * dimz, true, 2 * haloSize); break;
			}
		}

		ScalarField3D(BackendType _hw, ScalarField3D *field) : 
			hw(_hw), dimx(field->dimx), dimy(field->dimy), dimz(field->dimz),
			dx(field->dx), dy(field->dy), dz(field->dz), haloSize(field->haloSize)
		{
			PARAplan *pplan = PARAplan::Instance();
			dimxOffset = pplan->getOffset1D();
			switch( hw )
			{
			case CPU: 
				u = new FTYPE[dimx * dimy * dimz + 2 * haloSize]; 
				switch( field->hw )
				{
				case CPU: memcpy(u + haloSize, field->getArray() + haloSize, dimx * dimy * dimz * sizeof(FTYPE)); break;
				case GPU: multiDevMemcpy<FTYPE>(u + haloSize, field->getMultiArray(), dimx * dimy * dimz, haloSize); break;
				}
				break;
			case GPU: 
				multiDevAlloc<FTYPE> (dd_u, dimx * dimy * dimz, true, 2 * haloSize); 
				switch( field->hw )
				{
				case CPU: multiDevMemcpy<FTYPE>(dd_u, field->getArray() + field->haloSize, dimx * dimy * dimz, haloSize);  break;
				case GPU: multiDevMemcpy<FTYPE>(dd_u, field->getMultiArray(), dimx * dimy * dimz, haloSize); break;
				}
				break;
			}
		}

		~ScalarField3D()
		{
			switch( hw )
			{
			case CPU: delete [] u; break;
			case GPU: multiDevFree<FTYPE>(dd_u); break;
			}
		}

		void CopyFieldTo(Grid3D *grid, ScalarField3D *dest, NodeType type)
		{
			switch( hw )
			{
			case CPU:
				{
					for (int i = 0; i < dimx; i++)
						for (int j = 0; j < dimy; j++)
							for (int k = 0; k < dimz; k++)
								if (grid->GetType(i + dimxOffset, j, k) == type)
									dest->elem(i, j, k) = elem(i, j, k);
					break;
				}
			case GPU: 
				{
					CopyFieldTo_GPU(dimx, dimy, dimz, dd_u, dest->getMultiArray(), grid->GetTypesGPU(), type, haloSize); 
					break;
				}
			}
		}

	void MergeFieldTo(Node *nodes, ScalarField3D *dest, NodeType type)
	{
		switch( hw )
		{
		case CPU:
			{
				#pragma omp parallel default(none) firstprivate(type) shared(nodes, dest)
				{
					#pragma omp for
					for (int i = 0; i < dimx; i++)
						for (int j = 0; j < dimy; j++)
							for (int k = 0; k < dimz; k++)
							{
								int id = (i + dimxOffset) * dimy * dimz + j * dimz + k;
								if (nodes[id].type == type)
									dest->elem(i, j, k) = (dest->elem(i, j, k) + elem(i, j, k)) / 2;
							}
				}
				break;
			}
		}
	}

		void MergeFieldTo(NodeType **nodes, ScalarField3D *dest, NodeType type)
		{
			switch( hw )
			{
			case GPU:
				{
					MergeFieldTo_GPU(dimx, dimy, dimz, dd_u, dest->getMultiArray(), nodes, type, haloSize);
					break;
				}
			}
		}

		void Smooth(Node *nodes, ScalarField3D *dest, NodeType type)
		{
			switch( hw )
			{
			case CPU:
				{
					#pragma omp parallel default(none) firstprivate(type) shared(nodes, dest)
					{
						#pragma omp for
						for (int i = 0; i < dimx; i++)
							for (int j = 0; j < dimy; j++)
								for (int k = 0; k < dimz; k++)
								{
									int id = i * dimy * dimz + j * dimz + k;
									if (nodes[id].type == type)
										dest->elem(i, j, k) = (elem(i, j, k) + elem(i+1, j, k) + elem(i-1, j, k) + 
															   elem(i, j-1, k) + elem(i, j+1, k) + 
															   elem(i, j, k-1) + elem(i, j, k+1)) / 7;
								}
					}
					break;
				}
			}
		}

		void Smooth(Node **nodes, ScalarField3D *dest, NodeType type)
		{
			switch( hw )
			{
			case GPU:
				{
					throw logic_error("not implemented yet!\n");
				}
			}
		}

		void Transpose(ScalarField3D *dest)
		{
			switch( hw )
			{
			case CPU:
				{
					throw logic_error("Transpose is not implemented for CPU!\n");
				}
			case GPU:
				{
					Transpose_GPU(dimx, dimy, dimz, getMultiArray(), dest->getMultiArray(), haloSize);
					break;	
				}
			}		
		}

		void DumpToFile(char *filename, int x = -1)
		{
			PARAplan *pplan = PARAplan::Instance();
			if (pplan->size() > 1)
				throw runtime_error("DumpToFile() is not supported in MPI version");
			ScalarField3D host_copy(CPU, this);

			int x_start, x_end;
			if( x == -1 ) { x_start = 0; x_end = dimx; }
				else { x_start = x; x_end = x+1; }
			
			FILE *file = NULL;
			fopen_s(&file, filename, "w");
			for( int i = x_start; i < x_end; i++ )
			{
				fprintf(file, "x = %i\n", i);
				for( int j = 0; j < dimy; j++ )
				{
					for( int k = 0; k < dimz; k++ )
						fprintf(file, "%.3f ", host_copy.elem(i, j, k));
					fprintf(file, "\n");
				}
			}
			fclose(file);
		}

	private:
		int haloSize;
		int dimxOffset;
		FTYPE *u;		// store on CPU
		FTYPE **dd_u; // or multiGPU according to the hw flag
	};


	struct TimeLayer3D
	{
		BackendType hw;
		int dimx, dimy, dimz;
		FTYPE dx, dy, dz;

		ScalarField3D *U, *V, *W, *T;

		int haloSize;

		void syncHalos(FTYPE *mpi_buf = NULL)
		{
			U->syncHalos(666, 667, mpi_buf);
			V->syncHalos(666, 667, mpi_buf);
			W->syncHalos(666, 667, mpi_buf);
			T->syncHalos(666, 667, mpi_buf);
		}
		
		inline FTYPE DissFuncX(int i, int j, int k)
		{
			FTYPE u_x = U->d_x(i, j, k);
			FTYPE v_x = V->d_x(i, j, k);
			FTYPE w_x = W->d_x(i, j, k);

			FTYPE u_y = U->d_y(i, j, k);
			FTYPE u_z = U->d_z(i, j, k);

			return 2 * u_x * u_x + v_x * v_x + w_x * w_x + v_x * u_y + w_x * u_z;
		}

		inline FTYPE DissFuncY(int i, int j, int k)
		{
			FTYPE u_y = U->d_y(i, j, k);
			FTYPE v_y = V->d_y(i, j, k);
			FTYPE w_y = W->d_y(i, j, k);

			FTYPE v_x = V->d_x(i, j, k);
			FTYPE v_z = V->d_z(i, j, k);

			return u_y * u_y + 2 * v_y * v_y + w_y * w_y + u_y * v_x + w_y * v_z;
		}

		inline FTYPE DissFuncZ(int i, int j, int k)
		{
			FTYPE u_z = U->d_z(i, j, k);
			FTYPE v_z = V->d_z(i, j, k);
			FTYPE w_z = W->d_z(i, j, k);

			FTYPE w_x = W->d_x(i, j, k);
			FTYPE w_y = W->d_y(i, j, k);

			return u_z * u_z + v_z * v_z + 2 * w_z * w_z + u_z * w_x + v_z * w_y;
		}

		inline FTYPE DissFunc(int i, int j, int k)
		{
			return DissFuncX(i, j, k) + DissFuncY(i, j, k) + DissFunc(i, j, k);
		}

		double EvalDivError(Grid3D *grid)
		{
			ScalarField3D *U_cpu = new ScalarField3D(CPU, U);
			ScalarField3D *V_cpu = new ScalarField3D(CPU, V);
			ScalarField3D *W_cpu = new ScalarField3D(CPU, W);

			U_cpu->syncHalos(666, 667);
			V_cpu->syncHalos(668, 669);
			W_cpu->syncHalos(670, 671);

			PARAplan* pplan = PARAplan::Instance();
			int ndimx = (pplan->rank() == pplan->size()-1)? dimx-1:dimx;
			double err = 0.0;
			int count = 0;
			for (int i = 0; i < ndimx; i++)
				for (int j = 0; j < dimy-1; j++)
					for (int k = 0; k < dimz-1; k++)
						if (grid->GetType(i + dimxOffset, j, k) == NODE_IN)
						{
							double err_x = (U_cpu->elem(i, j, k) + U_cpu->elem(i, j-1, k) + U_cpu->elem(i, j-1, k-1) + U_cpu->elem(i, j, k-1) -
								U_cpu->elem(i-1, j, k) - U_cpu->elem(i-1, j-1, k) - U_cpu->elem(i-1, j-1, k-1) - U_cpu->elem(i-1, j, k-1)) * dz * dy / 4.0;

							double err_y = (V_cpu->elem(i, j, k) + V_cpu->elem(i-1, j, k) + V_cpu->elem(i-1, j, k-1) + V_cpu->elem(i, j, k-1) -
								V_cpu->elem(i, j-1, k) - V_cpu->elem(i-1, j-1, k) - V_cpu->elem(i-1, j-1, k-1) - V_cpu->elem(i, j-1, k-1)) * dx * dz / 4.0;

							double err_z = (W_cpu->elem(i, j, k) + W_cpu->elem(i, j-1, k) + W_cpu->elem(i-1, j-1, k) + W_cpu->elem(i-1, j, k) -
								W_cpu->elem(i, j, k-1) - W_cpu->elem(i, j-1, k-1) - W_cpu->elem(i-1, j-1, k-1) - W_cpu->elem(i-1, j, k-1)) * dx * dy / 4.0;

							err += abs(err_x + err_y + err_z);
							count++;
						}
			delete U_cpu;
			delete V_cpu;
			delete W_cpu;

#ifdef __PARA
			double err_total;
			int count_total;
			MPI_Reduce(&err, &err_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Bcast(&err_total, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Reduce(&count, &count_total, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Bcast(&count_total, 1, MPI_INT, 0, MPI_COMM_WORLD);
			return err_total / count_total;
#else
			return err / count;
#endif
		}

		void Smooth(Grid3D *grid, TimeLayer3D *dest, NodeType type)
		{
			Node *cpu_nodes = grid->GetNodesCPU();
			Node **dev_nodes = NULL;// = grid->GetNodesGPU(false);
			switch (hw)
			{
			case CPU:
				U->Smooth( cpu_nodes, dest->U, type );
				V->Smooth( cpu_nodes, dest->V, type );
				W->Smooth( cpu_nodes, dest->W, type );
				T->Smooth( cpu_nodes, dest->T, type );				
				break;
			case GPU:
				U->Smooth( dev_nodes, dest->U, type );
				V->Smooth( dev_nodes, dest->V, type );
				W->Smooth( dev_nodes, dest->W, type );
				T->Smooth( dev_nodes, dest->T, type );
				break;
			}
		}

		void MergeLayerTo(Grid3D *grid, TimeLayer3D *dest, NodeType type, bool transposed = false)
		{
			Node *cpu_nodes = grid->GetNodesCPU();
			NodeType **dev_nodes = grid->GetTypesGPU(transposed);
			//Node** dev_nodes = grid->GetNodesGPU();
			switch (hw)
			{
			case CPU:
				U->MergeFieldTo( cpu_nodes, dest->U, type );
				V->MergeFieldTo( cpu_nodes, dest->V, type );
				W->MergeFieldTo( cpu_nodes, dest->W, type );
				T->MergeFieldTo( cpu_nodes, dest->T, type );
				break;
			case GPU:
				U->MergeFieldTo( dev_nodes, dest->U, type );
				V->MergeFieldTo( dev_nodes, dest->V, type );
				W->MergeFieldTo( dev_nodes, dest->W, type );
				T->MergeFieldTo( dev_nodes, dest->T, type );
			}
		}

		void CopyLayerTo(TimeLayer3D *dest)
		{
			switch( hw )
			{
			case CPU:
				switch( dest->hw )
				{
				case CPU:
					memcpy(dest->U->getArray() + dest->haloSize, U->getArray() + haloSize, dimx * dimy * dimz * sizeof(FTYPE));
					memcpy(dest->V->getArray() + dest->haloSize, V->getArray() + haloSize, dimx * dimy * dimz * sizeof(FTYPE));
					memcpy(dest->W->getArray() + dest->haloSize, W->getArray() + haloSize, dimx * dimy * dimz * sizeof(FTYPE));
					memcpy(dest->T->getArray() + dest->haloSize, T->getArray() + haloSize, dimx * dimy * dimz * sizeof(FTYPE));
					return;
				case GPU:
					multiDevMemcpy<FTYPE>(dest->U->getMultiArray(), U->getArray() + haloSize, dimx * dimy * dimz, haloSize); 
					multiDevMemcpy<FTYPE>(dest->V->getMultiArray(), V->getArray() + haloSize, dimx * dimy * dimz, haloSize); 
					multiDevMemcpy<FTYPE>(dest->W->getMultiArray(), W->getArray() + haloSize, dimx * dimy * dimz, haloSize); 
					multiDevMemcpy<FTYPE>(dest->T->getMultiArray(), T->getArray() + haloSize, dimx * dimy * dimz, haloSize); 
					return;
				}
				break;
			case GPU:
				switch( dest->hw )
				{
				case CPU:
					multiDevMemcpy<FTYPE>(dest->U->getArray() + dest->haloSize, U->getMultiArray(), dimx * dimy * dimz, haloSize);
					multiDevMemcpy<FTYPE>(dest->V->getArray() + dest->haloSize, V->getMultiArray(), dimx * dimy * dimz, haloSize);
					multiDevMemcpy<FTYPE>(dest->W->getArray() + dest->haloSize, W->getMultiArray(), dimx * dimy * dimz, haloSize);
					multiDevMemcpy<FTYPE>(dest->T->getArray() + dest->haloSize, T->getMultiArray(), dimx * dimy * dimz, haloSize);
					return;
				case GPU:
					multiDevMemcpy<FTYPE>(dest->U->getMultiArray(), U->getMultiArray(), dimx * dimy * dimz, haloSize);
					multiDevMemcpy<FTYPE>(dest->V->getMultiArray(), V->getMultiArray(), dimx * dimy * dimz, haloSize);
					multiDevMemcpy<FTYPE>(dest->W->getMultiArray(), W->getMultiArray(), dimx * dimy * dimz, haloSize);
					multiDevMemcpy<FTYPE>(dest->T->getMultiArray(), T->getMultiArray(), dimx * dimy * dimz, haloSize);
					return;
				}
				break;
			}
		}

		void CopyLayerTo(Grid3D *grid, TimeLayer3D *dest, NodeType type)
		{			
			U->CopyFieldTo(grid, dest->U, type);
			V->CopyFieldTo(grid, dest->V, type);
			W->CopyFieldTo(grid, dest->W, type);
			T->CopyFieldTo(grid, dest->T, type);
		}

		void CopyFromGrid(Grid3D *grid)
		{
			switch( hw )
			{
			case CPU: 
				{
					for (int i = 0; i < dimx; i++)
						for (int j = 0; j < dimy; j++)
							for (int k = 0; k < dimz; k++)
							{
								Vec3D vel = grid->GetVel(i + dimxOffset, j, k);
								U->elem(i, j, k) = vel.x;
								V->elem(i, j, k) = vel.y;
								W->elem(i, j, k) = vel.z;
								T->elem(i, j, k) = (FTYPE)grid->GetT(i + dimxOffset, j, k);
							}
					break;
				}
				case GPU:
					{
						FTYPE *u_cpu = new FTYPE[dimx * dimy * dimz];
						FTYPE *v_cpu = new FTYPE[dimx * dimy * dimz];
						FTYPE *w_cpu = new FTYPE[dimx * dimy * dimz];
						FTYPE *T_cpu = new FTYPE[dimx * dimy * dimz];
						for (int i = 0; i < dimx; i++)
							for (int j = 0; j < dimy; j++)
								for (int k = 0; k < dimz; k++)
								{
									int id = i * dimy * dimz + j * dimz + k;
									Vec3D vel = grid->GetVel(i + dimxOffset, j, k);
									u_cpu[id] = vel.x;
									v_cpu[id] = vel.y;
									w_cpu[id] = vel.z;
									T_cpu[id] = grid->GetT(i + dimxOffset, j, k);
								}
						multiDevMemcpy<FTYPE>(U->getMultiArray(), u_cpu, dimx * dimy * dimz, haloSize);
						multiDevMemcpy<FTYPE>(V->getMultiArray(), v_cpu, dimx * dimy * dimz, haloSize);
						multiDevMemcpy<FTYPE>(W->getMultiArray(), w_cpu, dimx * dimy * dimz, haloSize);
						multiDevMemcpy<FTYPE>(T->getMultiArray(), T_cpu, dimx * dimy * dimz, haloSize); 

						delete [] u_cpu;
						delete [] v_cpu;
						delete [] w_cpu;
						delete [] T_cpu;
						break;
					}
			}
		}

		void CopyToGrid(Grid3D *grid)
		{
			switch( hw )
			{
			case CPU:
				{
					for (int i = 0; i < dimx; i++)
						for (int j = 0; j < dimy; j++)
							for (int k = 0; k < dimz; k++)
								grid->SetNodeVel(i + dimxOffset, j, k, Vec3D(U->elem(i, j, k), V->elem(i, j, k), W->elem(i, j , k)));
					break;
				}
				case GPU:
				{
					FTYPE *u_cpu = new FTYPE[dimx * dimy * dimz];
					FTYPE *v_cpu = new FTYPE[dimx * dimy * dimz];
					FTYPE *w_cpu = new FTYPE[dimx * dimy * dimz];
					multiDevMemcpy<FTYPE>(u_cpu, U->getMultiArray(), dimx * dimy * dimz, haloSize);
					multiDevMemcpy<FTYPE>(v_cpu, V->getMultiArray(), dimx * dimy * dimz, haloSize);
					multiDevMemcpy<FTYPE>(w_cpu, W->getMultiArray(), dimx * dimy * dimz, haloSize);
					for (int i = 0; i < dimx; i++)
						for (int j = 0; j < dimy; j++)
							for (int k = 0; k < dimz; k++)
								grid->SetNodeVel(i + dimxOffset, j, k, Vec3D(u_cpu[i * dimy * dimz + j * dimz + k],
																v_cpu[i * dimy * dimz + j * dimz + k],
																w_cpu[i * dimy * dimz + j * dimz + k]));

					delete [] u_cpu;
					delete [] v_cpu;
					delete [] w_cpu;

					break;
				}
			}
		}

		void FilterToArrays(Vec3D *outV, double *outT, int _outdimx, int outdimy, int outdimz)
		{
			if (_outdimx == 0) _outdimx = dimx;
			if (outdimy == 0) outdimy = dimy;
			if (outdimz == 0) outdimz = dimz;

			PARAplan *pplan = PARAplan::Instance();
			int outdimx, offset;
			pplan->get1D(outdimx, offset, _outdimx);

			int size = outdimy * outdimz;
			size *= (pplan->rank()==0)? _outdimx:outdimx;

			FTYPE *outVx = new FTYPE[size];
			FTYPE *outVy = new FTYPE[size];
			FTYPE *outVz = new FTYPE[size];

			switch( hw )
			{
			case CPU:
				{
					if (pplan->size() > 1)
						throw logic_error("MPI version of FilterToArrays is not implemented for 'CPU' backend");
					for (int i = 0; i < outdimx; i++)
						for (int j = 0; j < outdimy; j++)
							for (int k = 0; k < outdimz; k++)
							{	
								int x = (i * dimx / outdimx);
								int y = (j * dimy / outdimy);
								int z = (k * dimz / outdimz);
								int ind = i * outdimy * outdimz + j * outdimz + k;
								outVx[ind] = U->elem(x, y, z); 
								outVy[ind] = V->elem(x, y, z);
								outVz[ind] = W->elem(x, y, z);
								outT[ind] = T->elem(x, y, z);
							}
					break;
				}
			case GPU:
				{
					FTYPE *u_cpu = new FTYPE[dimx * dimy * dimz];
					FTYPE *v_cpu = new FTYPE[dimx * dimy * dimz];
					FTYPE *w_cpu = new FTYPE[dimx * dimy * dimz];
					FTYPE *T_cpu = new FTYPE[dimx * dimy * dimz];

					multiDevMemcpy<FTYPE>(u_cpu, U->getMultiArray(), dimx * dimy * dimz, haloSize);
					multiDevMemcpy<FTYPE>(v_cpu, V->getMultiArray(), dimx * dimy * dimz, haloSize);
					multiDevMemcpy<FTYPE>(w_cpu, W->getMultiArray(), dimx * dimy * dimz, haloSize);
					multiDevMemcpy<FTYPE>(T_cpu, T->getMultiArray(), dimx * dimy * dimz, haloSize);

					for (int i = 0; i < outdimx; i++)
						for (int j = 0; j < outdimy; j++)
							for (int k = 0; k < outdimz; k++)
							{	
								int x = (i * dimx / outdimx);
								int y = (j * dimy / outdimy);
								int z = (k * dimz / outdimz);
								int ind = i * outdimy * outdimz + j * outdimz + k;
								outVx[ind] = u_cpu[x * dimy * dimz + y * dimz + z]; 
								outVy[ind] = v_cpu[x * dimy * dimz + y * dimz + z]; 
								outVz[ind] = w_cpu[x * dimy * dimz + y * dimz + z]; 
								outT[ind] = T_cpu[x * dimy * dimz + y * dimz + z]; 
							}

					delete [] u_cpu;
					delete [] v_cpu;
					delete [] w_cpu;
					delete [] T_cpu;
					break;
				}
			}

#ifdef __PARA
			MPI_Status status;
			if (pplan->rank() > 0)
			{
				MPI_Send(&size, 1, MPI_INT, 0, 60, MPI_COMM_WORLD);
				MPI_Send(outVx, size, mpi_typeof(outVx), 0, 61, MPI_COMM_WORLD);
				MPI_Send(outVy, size, mpi_typeof(outVy), 0, 62, MPI_COMM_WORLD);
				MPI_Send(outVz, size, mpi_typeof(outVz), 0, 63, MPI_COMM_WORLD);
				MPI_Send(outT,  size, mpi_typeof(outT),  0, 64, MPI_COMM_WORLD);
			}			
			if (pplan->rank() == 0)
			{
				size = outdimx * outdimy * outdimz; offset = 0;
				for(int irank = 1; irank < pplan->size(); irank++)
				{
					offset += size;
					MPI_Recv(&size, 1, MPI_INT, irank, 60, MPI_COMM_WORLD, &status);
					fflush(stdout);
					MPI_Recv(outVx + offset, size, mpi_typeof(outVx), irank, 61, MPI_COMM_WORLD, &status);
					MPI_Recv(outVy + offset, size, mpi_typeof(outVy), irank, 62, MPI_COMM_WORLD, &status);
					MPI_Recv(outVz + offset, size, mpi_typeof(outVz), irank, 63, MPI_COMM_WORLD, &status);
					MPI_Recv(outT + offset, size,  mpi_typeof(outT),  irank, 64, MPI_COMM_WORLD, &status);
				}
			}
#endif
			if (pplan->rank() == 0)
				for (int i = 0; i < _outdimx * outdimy * outdimz; i++)
				{
					outV[i].x = outVx[i]; outV[i].y = outVy[i]; outV[i].z = outVz[i];
				}
			delete [] outVx;
			delete [] outVy;
			delete [] outVz;
		}

		void CopyFromGrid(Grid3D *grid, NodeType target)
		{
			switch( hw )
			{
			case CPU:
				{
					for (int i = 0; i < dimx; i++)
						for (int j = 0; j < dimy; j++)
							for (int k = 0; k < dimz; k++)
								if (grid->GetType(i + dimxOffset, j, k) == target)
								{
									Vec3D velocity = grid->GetVel(i + dimxOffset, j, k);
									U->elem(i, j, k) = (FTYPE)velocity.x;
									V->elem(i, j, k) = (FTYPE)velocity.y;
									W->elem(i, j, k) = (FTYPE)velocity.z;
									T->elem(i, j, k) = (FTYPE)grid->GetT(i + dimxOffset, j, k);
								}
					break;
				}
			case GPU:
				{
					//CopyFromGrid_GPU(dimx, dimy, dimz, U->getMultiArray(), V->getMultiArray(), W->getMultiArray(), T->getMultiArray(), grid->GetNodesGPU(), target, haloSize); 
					break;
				}
			}
		}

		void CopyGridBoundary(Grid3D *grid)
		{
			switch (hw)
			{
			case CPU:
				CopyFromGrid(grid, NODE_BOUND);
				CopyFromGrid(grid, NODE_VALVE);
				break;
			}			
		}

		void CopyGridBoundary(DirType dir, int *num_seg, Segment3D **segs,  NodesBoundary3D **nodes)
		{
			switch (hw)
			{
			case GPU:
				CopyGridBoundary_GPU(dir, dimx, dimy, dimz, U->getMultiArray(), V->getMultiArray(), W->getMultiArray(), T->getMultiArray(), num_seg, segs, nodes, haloSize);
				break;
			}			
		}

		void Clear(Grid3D *grid, NodeType target, FTYPE const_u, FTYPE const_v, FTYPE const_w, FTYPE const_T)
		{
			switch( hw )
			{
			case CPU:
				{
					for (int i = 0; i < dimx; i++)
						for (int j = 0; j < dimy; j++)
							for (int k = 0; k < dimz; k++)
								if (grid->GetType(i + dimxOffset, j, k) == target)
								{
									U->elem(i, j, k) = const_u;
									V->elem(i, j, k) = const_v;
									W->elem(i, j, k) = const_w;
									T->elem(i, j, k) = const_T;
								}
					break;
				}
			case GPU:
				{
					Clear_GPU(dimx, dimy, dimz, U->getMultiArray(), V->getMultiArray(), W->getMultiArray(), T->getMultiArray(), grid->GetTypesGPU(), target, const_u, const_v, const_w, const_T, haloSize);
					break;
				}
			}
		}

		void Transpose(TimeLayer3D *dest)
		{
			U->Transpose(dest->U);
			V->Transpose(dest->V);
			W->Transpose(dest->W);
			T->Transpose(dest->T);
		}

		void PrintToFile(FILE *file, const char *desc, FTYPE *u, int dimx, int dimy, int dimz)
		{
			PARAplan *pplan = PARAplan::Instance();
			if (pplan->rank() > 0)
				return;
			fprintf(file, "Array %s[%i,%i,%i]:\n", desc, dimx, dimy, dimz);
			for (int i = 0; i < dimx; i++)
			{			
				for (int j = 0; j < dimy; j++)
				{
					for (int k = 0; k < dimz; k++)
					{
						int id = i * dimy * dimz + j * dimz + k;
						fprintf(file, "%.8f ", u[id]);
					}
					fprintf(file, "\n");
				}
				fprintf(file, "\n");
			}
		}

		void PrintToFile(const char *filename)
		{
			PARAplan *pplan = PARAplan::Instance();
			if (pplan->rank() > 0)
				return;
			FILE *file = NULL;
			fopen_s(&file, filename, "w");

			if( hw == GPU )
			{
				FTYPE *u_cpu = new FTYPE[dimx * dimy * dimz];
				FTYPE *v_cpu = new FTYPE[dimx * dimy * dimz];
				FTYPE *w_cpu = new FTYPE[dimx * dimy * dimz];
				FTYPE *T_cpu = new FTYPE[dimx * dimy * dimz];

				multiDevMemcpy<FTYPE>(u_cpu, U->getMultiArray(), dimx * dimy * dimz, haloSize);
				multiDevMemcpy<FTYPE>(v_cpu, V->getMultiArray(), dimx * dimy * dimz, haloSize);
				multiDevMemcpy<FTYPE>(w_cpu, W->getMultiArray(), dimx * dimy * dimz, haloSize);
				multiDevMemcpy<FTYPE>(T_cpu, T->getMultiArray(), dimx * dimy * dimz, haloSize);

				PrintToFile(file, "u", u_cpu, dimx, dimy, dimz);
				PrintToFile(file, "v", v_cpu, dimx, dimy, dimz);
				PrintToFile(file, "w", w_cpu, dimx, dimy, dimz);
				PrintToFile(file, "T", T_cpu, dimx, dimy, dimz);
				
				delete [] u_cpu;
				delete [] v_cpu;
				delete [] w_cpu;
				delete [] T_cpu;
			}
			else
				printf("Print to file is not implemented for CPU! :)\n");

			fclose(file);
		}
		
		TimeLayer3D(BackendType _hw, int _dimx, int _dimy, int _dimz, FTYPE _dx, FTYPE _dy, FTYPE _dz, int _haloSize = 0) : 
			hw(_hw), dimx(_dimx), dimy(_dimy), dimz(_dimz),
			dx(_dx), dy(_dy), dz(_dz), haloSize(_haloSize)
		{
			PARAplan *pplan = PARAplan::Instance();
			dimxOffset = pplan->getOffset1D();
			U = new ScalarField3D(hw, dimx, dimy, dimz, dx, dy, dz, haloSize);
			V = new ScalarField3D(hw, dimx, dimy, dimz, dx, dy, dz, haloSize);
			W = new ScalarField3D(hw, dimx, dimy, dimz, dx, dy, dz, haloSize);
			T = new ScalarField3D(hw, dimx, dimy, dimz, dx, dy, dz, haloSize);
		}

		TimeLayer3D(BackendType _hw, Grid3D *grid, int _haloSize = 0) : 
			hw(_hw), dimx(grid->dimx), dimy(grid->dimy), dimz(grid->dimz),
			dx((FTYPE)grid->dx), dy((FTYPE)grid->dy), dz((FTYPE)grid->dz), haloSize(_haloSize)
		{
			PARAplan *pplan = PARAplan::Instance();
			dimx = pplan->getLength1D();
			dimxOffset = pplan->getOffset1D();
			U = new ScalarField3D(hw, dimx, dimy, dimz, dx, dy, dz, haloSize);
			V = new ScalarField3D(hw, dimx, dimy, dimz, dx, dy, dz, haloSize);
			W = new ScalarField3D(hw, dimx, dimy, dimz, dx, dy, dz, haloSize);
			T = new ScalarField3D(hw, dimx, dimy, dimz, dx, dy, dz, haloSize);

			CopyFromGrid(grid);
		}
		
		~TimeLayer3D()
		{
			delete U;
			delete V;
			delete W;
			delete T;
		}
	private:
		int dimxOffset;
	};

	struct TimeLayer3D_GPU
	{
		FTYPE *u, *v, *w, *T;
		FTYPE **dd_u, **dd_v, **dd_w, **dd_T;
		int dimx, dimy, dimz;
		int haloSize;

		inline __device__ FTYPE& elem(FTYPE *arr, int i, int j, int k)
		{
			return arr[haloSize + i * dimy * dimz + j * dimz + k];
		}

		inline __device__ FTYPE d_x(FTYPE *arr, FTYPE dx, int i, int j, int k)	{ return (elem(arr, i+1, j, k) - elem(arr, i-1, j, k)) / (2 * dx); }
		inline __device__ FTYPE d_y(FTYPE *arr, FTYPE dy, int i, int j, int k)	{ return (elem(arr, i, j+1, k) - elem(arr, i, j-1, k)) / (2 * dy); }
		inline __device__ FTYPE d_z(FTYPE *arr, FTYPE dz, int i, int j, int k)	{ return (elem(arr, i, j, k+1) - elem(arr, i, j, k-1)) / (2 * dz); }
		inline __device__ FTYPE d_xx(FTYPE *arr, FTYPE dx, int i, int j, int k)	{ return (elem(arr, i+1, j, k) - 2 * elem(arr, i, j, k) + elem(arr, i-1, j, k)) / (dx * dx); }
		inline __device__ FTYPE d_yy(FTYPE *arr, FTYPE dy, int i, int j, int k)	{ return (elem(arr, i, j+1, k) - 2 * elem(arr, i, j, k) + elem(arr, i, j-1, k)) / (dy * dy); }
		inline __device__ FTYPE d_zz(FTYPE *arr, FTYPE dz, int i, int j, int k)	{ return (elem(arr, i, j, k+1) - 2 * elem(arr, i, j, k) + elem(arr, i, j, k-1)) / (dz * dz); }

		inline __device__ FTYPE DissFuncX(FTYPE dx, FTYPE dy, FTYPE dz, int i, int j, int k)
		{
			FTYPE u_x = d_x(u, dx, i, j, k);
			FTYPE v_x = d_x(v, dx, i, j, k);
			FTYPE w_x = d_x(w, dx, i, j, k);

			FTYPE u_y = d_y(u, dy, i, j, k);
			FTYPE u_z = d_z(u, dz, i, j, k);

			return 2 * u_x * u_x + v_x * v_x + w_x * w_x + v_x * u_y + w_x * u_z;
		}

		inline __device__ FTYPE DissFuncY(FTYPE dx, FTYPE dy, FTYPE dz, int i, int j, int k)
		{
			FTYPE u_y = d_y(u, dy, i, j, k);
			FTYPE v_y = d_y(v, dy, i, j, k);
			FTYPE w_y = d_y(w, dy, i, j, k);

			FTYPE v_x = d_x(v, dx, i, j, k);
			FTYPE v_z = d_z(v, dz, i, j, k);

			return u_y * u_y + 2 * v_y * v_y + w_y * w_y + u_y * v_x + w_y * v_z;
		}

		inline __device__ FTYPE DissFuncZ(FTYPE dx, FTYPE dy, FTYPE dz, int i, int j, int k)
		{
			FTYPE u_z = d_z(u, dz, i, j, k);
			FTYPE v_z = d_z(v, dz, i, j, k);
			FTYPE w_z = d_z(w, dz, i, j, k);

			FTYPE w_x = d_x(w, dx, i, j, k);
			FTYPE w_y = d_y(w, dy, i, j, k);

			return u_z * u_z + v_z * v_z + 2 * w_z * w_z + u_z * w_x + v_z * w_y;
		}

		inline __device__ FTYPE DissFuncZ_as_Y(FTYPE dx, FTYPE dy, FTYPE dz, int i, int j, int k)
		{
			FTYPE u_z = d_y(u, dy, i, j, k);
			FTYPE v_z = d_y(v, dy, i, j, k);
			FTYPE w_z = d_y(w, dy, i, j, k);

			FTYPE w_x = d_x(w, dx, i, j, k);
			FTYPE w_y = d_z(w, dz, i, j, k);

			return u_z * u_z + v_z * v_z + 2 * w_z * w_z + u_z * w_x + v_z * w_y;
		}

	  void SetDevice(int i)
		{
			u = dd_u[i];
			v = dd_v[i];
			w = dd_w[i];
			T = dd_T[i];
		}

		TimeLayer3D_GPU(TimeLayer3D *layer)
		{
			u = layer->U->getArray();
			v = layer->V->getArray();
			w = layer->W->getArray();
			T = layer->T->getArray();

			dd_u = layer->U->getMultiArray();
			dd_v = layer->V->getMultiArray();
			dd_w = layer->W->getMultiArray();
			dd_T = layer->T->getMultiArray();

			dimx = layer->dimx;
			dimy = layer->dimy;
			dimz = layer->dimz;

			haloSize = layer->haloSize;
		}
	};
}
