#include "PARAplan.h"

namespace Common
{
#ifdef __PARA
	void mpiSafeCall(int status, char* message)
	{
		if (status != MPI_SUCCESS )
			throw std::runtime_error(message);		 
	}

	MPI_Datatype mpi_typeof(int*) { return MPI_INT; }
	MPI_Datatype mpi_typeof(float*) { return MPI_FLOAT; }
	MPI_Datatype mpi_typeof(double*) { return MPI_DOUBLE; }
#endif

	bool PARAplan::isInstance = false;
	PARAplan* PARAplan::self = NULL;

	PARAplan::PARAplan()
	{
		nNodes = 1; data1D = 0; iRank = 0;
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
#ifdef __PARA
			MPI_Reduce(&nGPU, &nGPUTotal, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Bcast(&nGPUTotal, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
#else
			nGPUTotal = nGPU;
#endif
		}
	}

	void PARAplan::getEven1D(int& split, int& offset, int num_elems_1D)
	{
		split = num_elems_1D / nNodes;
		offset = iRank*split;
		if (iRank < num_elems_1D % nNodes)
		{
			split++;
			offset += iRank;
		}
		else
			offset += num_elems_1D % nNodes;
	}

	void PARAplan::splitEven1D(int num_elems_1D)
	{
		getEven1D(data1D, offset1D, num_elems_1D);
		printf("PARAplan::splitEven1D: node %d: data1D = %d, offset1D = %d\n", iRank, data1D, offset1D);
		if (hw == GPU)
			pGPUplan->splitEven1D(data1D);
		fflush(stdout);
	}

	void PARAplan::setGPUnum(int num)
	{
		pGPUplan->setGPUnum(num);
		nGPU = pGPUplan->size();
#ifdef __PARA
			MPI_Reduce(&nGPU, &nGPUTotal, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Bcast(&nGPUTotal, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
#else
			nGPUTotal = nGPU;
#endif
	}

	int PARAplan::size()const
	{
		return nNodes;
	}

	int PARAplan::rank()const
	{
		return iRank;
	}

	int PARAplan::gpuNum()const
	{
		return nGPU;
	}

	int PARAplan::gpuTotal()const
	{
		return nGPUTotal;
	}

	PARAplan::~PARAplan()
	{
		delete pGPUplan;
#ifdef __PARA
		MPI_Finalize();
#endif
	}

}
