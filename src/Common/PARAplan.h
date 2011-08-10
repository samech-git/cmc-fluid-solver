#pragma once

#define __PARA

#ifdef __PARA
#include "mpi.h"
#endif

#include "GPUplan.h"
#include "Geometry.h" //for BackendType

namespace Common
{
	struct PARAplan
	{	// singleton class. should use "delete" explicitly for the destructor to be called
		static PARAplan* Instance();
		int getLength1D()const{return data1D;}
		int getOffset1D()const{return offset1D;}
		void init(BackendType backend);
		void setGPUnum(int num);
		int size()const;
		int rank()const;
		int gpuNum()const;
		int gpuTotal()const;
		void getEven1D(int& split, int& offset, int num_elems_1D);
		void splitEven1D(int num_elems_1D);
		~PARAplan();
		private:
			static bool isInstance;
			static PARAplan *self;
			PARAplan();
			int offset1D;
			int  nNodes;
			int data1D;
			int iRank;
			int nGPU;
			int nGPUTotal;
			BackendType hw;
			GPUplan* pGPUplan;
	};

#ifdef __PARA
	extern void mpiSafeCall(int status, char* message);

	template <typename T> MPI_Datatype mpi_typeof(T*) { return -1; }
	extern MPI_Datatype mpi_typeof(int*);
	extern MPI_Datatype mpi_typeof(float*);
	extern MPI_Datatype mpi_typeof(double*);
#endif
}
