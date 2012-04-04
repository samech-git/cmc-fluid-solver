#include "PARAplan.h"

namespace Common
{
	bool PARAplan::isInstance = false;
	PARAplan *PARAplan::self = NULL;

	PARAplan::PARAplan()
	{
		nNodes = 1; data1D = 0; iRank = 0; offset1D = 0; data1DTotal = 0; nGPU = 0; nGPUTotal = 0;
	}

  void PARAplan::init(BackendType backend = GPU)
	{
		hw = backend;
#ifdef __PARA
		MPI_Comm_size(MPI_COMM_WORLD, &nNodes);
		MPI_Comm_rank(MPI_COMM_WORLD, &iRank);
#endif
		if (backend == GPU)
		{
			pGPUplan = GPUplan::Instance();
			pGPUplan->init();
			nGPU = pGPUplan->size();
//-------------For testing: 1 GPU per rank------
//		pGPUplan->setBegin(iRank*nGPU);
//----------------------------------------------
#ifdef __PARA
			MPI_Reduce(&nGPU, &nGPUTotal, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Bcast(&nGPUTotal, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
#else
			nGPUTotal = nGPU;
#endif
		}
	}

	PARAplan* PARAplan::Instance()
	{
		if(!isInstance)
		{
			self = new PARAplan();
			isInstance = true;
			return self;
		}
		else
			return self;
	}

	void PARAplan::get1D(int& split, int& offset, int num_elems_1D)
	{
		double alpha = double(num_elems_1D) / data1DTotal;
		split = int (alpha * data1D);
		int num_elems_tmp = split;
#ifdef __PARA
		num_elems_tmp = 0;
		MPI_Reduce(&split, &num_elems_tmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Bcast(&num_elems_tmp, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
		if (iRank < num_elems_1D - num_elems_tmp)
			split++;
		offset = 0;
#ifdef __PARA
		int *splitting = new int[nNodes];
		MPI_Allgather(&split, 1, MPI_INT, splitting, 1, MPI_INT, MPI_COMM_WORLD);
		for (int i = 1; i <= iRank; i++)
			offset += splitting[i-1];
		delete [] splitting;
#endif
	}

	void PARAplan::splitEven1D(int num_elems_1D)
	{
		data1DTotal = num_elems_1D;
		data1D = num_elems_1D / nNodes;
		offset1D = iRank*data1D;
		if (iRank < num_elems_1D % nNodes)
		{
			data1D++;
			offset1D += iRank;
		}
		else
			offset1D += num_elems_1D % nNodes;
		printf("PARAplan::splitEven1D: node %d: data1D = %d, offset1D = %d\n", iRank, data1D, offset1D);
		if (hw == GPU)
			pGPUplan->splitEven1D(data1D);
		fflush(stdout);
	}

	void PARAplan::split1D(int *num_elems_1D)
	/*
		Num elements array of splittings per GPU, sorted
		by MPI node rank
	*/
	{
		if (hw == CPU) // No MPI 
		{
			data1D = data1DTotal = num_elems_1D[0];
			printf("PARAplan::split1D: node %d: data1D = %d, offset1D = %d\n", iRank, data1D, offset1D);
			return;
		}
		int *gpuPerNode = new int[nNodes];
#ifdef __PARA
		MPI_Allgather(&nGPU, 1, MPI_INT, gpuPerNode, 1, MPI_INT, MPI_COMM_WORLD);
#endif
		offset1D = 0;
		int i, j, gpu_cnt = 0;
		for (i = 0; i < iRank; i++)
			for(j = 0; j < gpuPerNode[i]; j++)
			{
				offset1D += num_elems_1D[gpu_cnt]; 
				gpu_cnt++;
			}

			data1D = 0;
			for (j = 0; j < nGPU; j++)
				data1D += num_elems_1D[gpu_cnt + j];
		
		data1DTotal = data1D;
#ifdef __PARA
		MPI_Reduce(&data1D, &data1DTotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Bcast(&data1DTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
		printf("PARAplan::split1D: node %d: data1D = %d, offset1D = %d\n", iRank, data1D, offset1D);
		pGPUplan->split1D(num_elems_1D + gpu_cnt);
		fflush(stdout);
	}

	void PARAplan::setGPUnum(int num)
	{
		pGPUplan->setGPUnum(num);
		nGPU = pGPUplan->size();
//-------------For testing: 1 GPU per rank------
			//pGPUplan->setBegin(iRank*nGPU);
//----------------------------------------------
#ifdef __PARA
			MPI_Reduce(&nGPU, &nGPUTotal, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Bcast(&nGPUTotal, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
#else
			nGPUTotal = nGPU;
#endif
	}

	PARAplan::~PARAplan()
	{
		if (hw == GPU)
			delete pGPUplan;
#ifdef __PARA
		MPI_Finalize();
#endif
	}

#ifdef __PARA
	void mpiSafeCall(int status, char *message)
	{
		if (status != MPI_SUCCESS )
			throw std::runtime_error(message);		 
	}

	MPI_Datatype mpi_typeof(int*) { return MPI_INT; }
	MPI_Datatype mpi_typeof(float*) { return MPI_FLOAT; }
	MPI_Datatype mpi_typeof(double*) { return MPI_DOUBLE; }
#endif
}
