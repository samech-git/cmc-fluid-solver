#pragma once

#include <cuda_runtime.h>
#include <cstdio>
#include <stdexcept>

#define MGPU

//for debugging purposes:
# define MGPU_EMU 1
#if MGPU_EMU
#define DEFAULT_DEVICE 0
#define cudaSetDevice(i) cudaSetDevice(DEFAULT_DEVICE)
#define GPU_NUM 16
#endif

namespace Common
{

	struct GPUNode
	{
		static const int nMaxEvents = 128; // max events per node
    cudaEvent_t event[nMaxEvents];
    cudaStream_t stream, stream2, stream3;

		void init();

		int getLength1D(){return data1D;}

		void setLength1D(int _data1D){ data1D = _data1D;}

		void destroy();

	private:
		int data1D; // length of grid data block per node along chosen direction
	};

	struct GPUplan
	{	// singleton class. should use "delete" explicitly for the destructor to be called

		void init();

		static GPUplan* Instance();
		GPUNode* node(int i);

		int getLength1D()const{return data1D;}
		int begin()const{ return ibegin;}
		int size()const{ return nGPU;};

		void setDevice(int iDev);
		void deviceSynchronize();
		void setBegin(int iDev);
		void setGPUnum(int num); // will use num < nMaxGPU GPUs
		void splitEven1D(int num_elems_1D);
		void split1D(int *num_elems_1D);
		inline int rescale1D(int inode, int num_elems_total) const 
			{ return plan[inode]->getLength1D() * (num_elems_total / data1D); } 		// returns num_elems_total per node according to GPUplan splitting along one dimension

		void destroy();
		~GPUplan();

		private:
			GPUNode** plan;
			static bool isInstance;
			static GPUplan *self;
			GPUplan();

			int nMaxGPU;
			int  nGPU;
			int ibegin;
			int data1D;
	};

	extern void gpuSafeCall(cudaError status, char* message, int id = -1, char *file = NULL, int linenum = -1);

	extern GPUplan *CheckGPUplan(int num_elems_total);

	template <typename T>
	void multiDevAlloc(T**& d_array, int num_elems_total, bool splitting = true, int num_elems_add = 0)
	/*
		Allocates num_elems_total on multiple devices,
		splitting them according to the GPUplan
		sliceWidth - dimension  of additional padding (2*haloSize) space along GPUplan direction (dimx) per GPU
	*/
	{
		if (!splitting)
		{ 
			GPUplan *pGPUplan = GPUplan::Instance();
			d_array =  new T*[pGPUplan->size()];
			for (int i = 0; i < pGPUplan->size(); i++)
			{
				pGPUplan->setDevice(i);
				gpuSafeCall( cudaMalloc(&(d_array[i]), sizeof(T) * (num_elems_total + num_elems_add)),"multiDevAlloc", i);
				//printf("multiDevAlloc: allocated %.2f MB on device  %d\n", sizeof(T) * (num_elems_total + num_elems_add)/(1024.f*1024.f), i);
			}
			return;
		}

		GPUplan *pGPUplan = CheckGPUplan(num_elems_total);
		d_array =  new T*[pGPUplan->size()];
		for (int i = 0; i < pGPUplan->size(); i++)
		{ 
			int num_elems = pGPUplan->rescale1D(i, num_elems_total);;
			pGPUplan->setDevice(i);
			gpuSafeCall( cudaMalloc(&(d_array[i]), sizeof(T) * (num_elems + num_elems_add)), "multiDevAlloc", i );
			//printf("multiDevAlloc: allocated %.2f MB on device  %d\n", sizeof(T) * (num_elems_total + num_elems_add)/(1024.f*1024.f), i);
		}
	}

	template <typename T>
	void multiDevFree(T** d_array)
	{
		GPUplan *pGPUplan = GPUplan::Instance();
		for (int i = 0; i < pGPUplan->size(); i++)
		{
			pGPUplan->setDevice(i);
			//printf("cudaMultiFree: deallocating device %d\n", i);
			gpuSafeCall(cudaFree(d_array[i]), "multiDevFree", i);
		}
		delete[] d_array;
	}

	template <typename T>
	void multiDevMemcpy(T** dev_dest_array, T* hst_src_array, int num_elems_total, int deviceOffset = 0)
	/*
		Copies num_elems_total of  hst_src_array to multiple devices
		and splits the data according to the GPUplan 
		deviceOffset - initial data position on each GPU node is shifted by a constant
	*/
	{
		GPUplan *pGPUplan = CheckGPUplan(num_elems_total);
		int offset = 0;
		for (int i = 0; i < pGPUplan->size(); i++)
		{
			int num_elems = pGPUplan->rescale1D(i, num_elems_total);
				pGPUplan->setDevice(i);
			gpuSafeCall( cudaMemcpyAsync(dev_dest_array[i] + deviceOffset, hst_src_array + offset, sizeof(T) * num_elems, cudaMemcpyHostToDevice, pGPUplan->node(i)->stream), "multiDevMemcpy hostToDevice", i);
			offset += num_elems;
		}
		pGPUplan->deviceSynchronize();
	}

	template <typename T>
	void multiDevMemcpy(T* hst_dest_array, T** dev_src_array, int num_elems_total, int deviceOffset = 0)
	/*
		Copies dev_src_array into "num_elems_total" of hst_dest_array 
		and assembles the data on host according to the GPUplan 
	*/
	{
		GPUplan *pGPUplan = CheckGPUplan(num_elems_total);
		int offset = 0;
		for (int i = 0; i < pGPUplan->size(); i++)
		{
			int num_elems = pGPUplan->rescale1D(i, num_elems_total);
			pGPUplan->setDevice(i);
			gpuSafeCall( cudaMemcpyAsync(hst_dest_array + offset, dev_src_array[i] + deviceOffset, sizeof(T) * num_elems, cudaMemcpyDeviceToHost, pGPUplan->node(i)->stream), "multiDevMemcpy deviceToHost", i );
			offset += num_elems;
		}
		pGPUplan->deviceSynchronize();
	}

	template <typename T>
	void multiDevMemcpy(T** dev_dest_array, T** dev_src_array, int num_elems_total, int deviceOffset = 0)
	/*
		Copies dev_src_array into "num_elems_total" of dev_dest_array, which 
		is split according to the GPUplan
	*/
	{
		GPUplan *pGPUplan = CheckGPUplan(num_elems_total);
		int offset = 0;
		for (int i = 0; i < pGPUplan->size(); i++)
		{
			int num_elems = pGPUplan->rescale1D(i, num_elems_total);
			pGPUplan->setDevice(i);
			gpuSafeCall( cudaMemcpyAsync(dev_dest_array[i] + deviceOffset, dev_src_array[i] + deviceOffset, sizeof(T) * num_elems, cudaMemcpyDeviceToDevice, pGPUplan->node(i)->stream), "multiDevMemcpy deviceToDevice", i );
			offset += num_elems;
		}
		pGPUplan->deviceSynchronize();
	}
}
