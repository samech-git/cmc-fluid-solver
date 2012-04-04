#pragma once

//#define __PARA

#ifdef __PARA
#include <mpi.h>
#endif

#include "GPUplan.h"
#include "Geometry.h" //for BackendType

namespace Common
{
	struct PARAplan
	{	//should use "delete" explicitly for the destructor to be called
		
		static PARAplan* Instance();

		void init(BackendType backend);

		int getLength1D()const{return data1D;}
		int getOffset1D()const{return offset1D;}
		int size()const{return nNodes;}
		int rank()const{return iRank;}
		int gpuNum()const{return nGPU;}
		int gpuTotal()const{return nGPUTotal;}
		void get1D(int& split, int& offset, int num_elems_1D);

		void splitEven1D(int num_elems_1D);
		void split1D(int *num_elems_1D);
		void setGPUnum(int num);

		~PARAplan();
		private:
			static PARAplan *self;
			static bool isInstance;
			BackendType hw;
			GPUplan* pGPUplan;
			PARAplan();
			int offset1D;
			int  nNodes;
			int data1D;
			int data1DTotal;
			int iRank;
			int nGPU;
			int nGPUTotal;
	};

#ifdef __PARA
	extern void mpiSafeCall(int status, char* message);

	template <typename T> MPI_Datatype mpi_typeof(T*) { return MPI_Datatype(-1); }
	extern MPI_Datatype mpi_typeof(int*);
	extern MPI_Datatype mpi_typeof(float*);
	extern MPI_Datatype mpi_typeof(double*);
#endif
}
