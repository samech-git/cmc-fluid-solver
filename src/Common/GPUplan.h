#pragma once

#include <cuda_runtime.h>
#include <stdio.h>
#include <stdexcept>
#include <cstring>
//#include "cmc_error.h" // for error
//#include <cmath>

#define MGPU

// for debugging purposes:
# define MGPU_EMU 0
#if MGPU_EMU
#define DEFAULT_DEVICE 0
	#define cudaSetDevice(i) cudaSetDevice(DEFAULT_DEVICE)
	#define GPU_NUM 16
#endif

namespace Common
{

	struct GPUNode
	{
		static const int nMaxEvents = 64; // max events per node
    cudaEvent_t event[nMaxEvents];

    //Stream for asynchronous command execution
    cudaStream_t stream, stream2, stream3;

		int getLength1D(){return data1D;}
		void setLength1D(int _data1D){ data1D = _data1D;}

		void init();
		void destroy();

	private:
		int data1D; // length of grid data block per node along chosen direction
	};

	struct GPUplan
	{	// singleton class. should use "delete" explicitly for the destructor to be called
		static GPUplan* Instance();
		GPUNode* node(int i);// {return this->operator[](i);}
		int getLength1D()const{return data1D;}
		// rescale1D() returns num_elems_total per node according to GPUplan splitting along one dimension
		inline int rescale1D(int inode, int num_elems_total) const { return plan[inode]->getLength1D() * (num_elems_total / data1D); }
		void init();
		void setGPUnum(int num); // will use num < nMaxGPU GPUs
		int size()const;
		void splitEven1D(int num_elems_1D);
		void destroy();
		~GPUplan();
		private:
			static bool isInstance;
			static GPUplan *self;
			GPUplan();
			int  nGPU;
			int nMaxGPU;
			GPUNode** plan;
			int data1D;
			void* hstRegPtr;
	};

	extern void gpuSafeCall(cudaError_t status, std::string message);

	extern GPUplan * CheckGPUplan(int num_elems_total);

	template <typename T>
	void multiDevAlloc(T**& d_array, int num_elems_total, bool splitting = true, int num_elems_add = 0)
	/*
		Allocates num_elems_total on multiple devices,
		splitting them according to the GPUplan
		sliceWidth - dimension  of additional padding (2*haloSize) space along GPUplan direction (dimx) per GPU
	*/
	{
		cudaError_t status = cudaSuccess;
		if (!splitting)
		{ 
			GPUplan *pGPUplan = GPUplan::Instance();
			d_array =  new T*[pGPUplan->size()];
			for (int i = 0; i < pGPUplan->size(); i++)
			{
				cudaSetDevice(i);
				//printf("multiDevAlloc: allocated %.2f MB on device  %d\n", sizeof(T) * (num_elems_total + num_elems_add)/(1024.f*1024.f), i);
				status = cudaMalloc(&(d_array[i]), sizeof(T) * (num_elems_total + num_elems_add));
				if (status != cudaSuccess )
					throw std::runtime_error("multiDevAlloc"); 	
			}
			return;
		}

		GPUplan *pGPUplan = CheckGPUplan(num_elems_total);
		d_array =  new T*[pGPUplan->size()];
		for (int i = 0; i < pGPUplan->size(); i++)
		{ 
			int num_elems = pGPUplan->rescale1D(i, num_elems_total);;
			cudaSetDevice(i);
			//printf("multiDevAlloc: allocated %.2f MB on device  %d\n", sizeof(T) * (num_elems_total + num_elems_add)/(1024.f*1024.f), i);
			status = cudaMalloc(&(d_array[i]), sizeof(T) * (num_elems + num_elems_add));
			if (status != cudaSuccess )
				throw std::runtime_error("multiDevAlloc"); 			
		}
	}

	template <typename T>
	void multiDevFree(T** d_array)
	{
		GPUplan *pGPUplan = GPUplan::Instance();
		for (int i = 0; i < pGPUplan->size(); i++)
		{
			cudaSetDevice((int)i);
			//printf("cudaMultiFree: deallocating device %d\n", i);
			cudaFree(d_array[i]);
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
			cudaSetDevice(i);
			//cudaMemcpy(dev_dest_array[i] + offset, hst_src_array + offset, sizeof(T) * num_elems, cudaMemcpyHostToDevice);
			cudaMemcpyAsync(dev_dest_array[i] + deviceOffset, hst_src_array + offset, sizeof(T) * num_elems, cudaMemcpyHostToDevice, pGPUplan->node(i)->stream);
			offset += num_elems;
		}
		cudaDeviceSynchronize();
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
			cudaSetDevice(i);
			//cudaMemcpy(hst_dest_array + offset, dev_src_array[i] + offset, sizeof(T) * num_elems, cudaMemcpyDeviceToHost);
			cudaMemcpyAsync(hst_dest_array + offset, dev_src_array[i] + deviceOffset, sizeof(T) * num_elems, cudaMemcpyDeviceToHost, pGPUplan->node(i)->stream);
			offset += num_elems;
		}
		cudaDeviceSynchronize();
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
			cudaSetDevice(i);
			//cudaMemcpy(dev_dest_array[i] + offset, dev_src_array[i] + offset, sizeof(T) * num_elems, cudaMemcpyDeviceToDevice);
			cudaMemcpyAsync(dev_dest_array[i] + deviceOffset, dev_src_array[i] + deviceOffset, sizeof(T) * num_elems, cudaMemcpyDeviceToDevice, pGPUplan->node(i)->stream);
			offset += num_elems;
		}
		cudaDeviceSynchronize();
	}
}
