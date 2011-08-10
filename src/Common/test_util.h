#pragma once

#include "../FluidSolver3D/Grid3D.h"
#include "GPUplan.h"
#include <cmath>

using namespace Common;

namespace TestUtil
{
	template<typename T>
	void fillRandom(T* arr, size_t num_elems)
	{
		for (size_t i = 0; i < num_elems; i++)
			arr[i] = T(10) + std::rand() % 512;
	}

	template<typename T>
	double sumEllements(T* arr, size_t num_elems)
	{
		double sum = 0.;
		for (size_t i = 0; i < num_elems; i++)
			sum += arr[i];
		return sum;
	};

	extern double sumEllements(FluidSolver3D::Node* arr, size_t num_elems);

	template<typename T>
	double sumEllementsMultiGPU(T** arr, size_t num_elems_total, size_t sliceOffset = 0)
	{
		double sum = 0.;
		T *arr_cpu = new T[num_elems_total];
		multiDevMemcpy<T>(arr_cpu, arr, num_elems_total, sliceOffset);
		sum = sumEllements<T>(arr_cpu,  num_elems_total);
		/* Debug:
		printf("sumEllementsMultiGPU: printing multiarray:\n");
		for (size_t i = 0; i < num_elems_total; i++)
			printf("%f  ", arr_cpu[i]);
		printf("\n");
		/**/
		delete [] arr_cpu;
		return sum;
	}

	extern double sumEllementsMultiGPU(FluidSolver3D::Node** arr, size_t num_elems_total);

	template <class T>
	void printEllementsMultiGPU(T** arr, int dimy, int dimz, int sliceOffset = 0, bool printWithOffset = false)
	{
		GPUplan* pGPUplan = GPUplan::Instance();
		int dimx = pGPUplan->getLength1D();
		dimx += (printWithOffset)?2*pGPUplan->size():0;
		//throw std::runtime_error("num_elems_total should be divisible! Terminating..."); 
		T *arr_cpu = new T[dimx*dimy*dimz];
		if (printWithOffset)
			multiDevMemcpy<T>(arr_cpu, arr, dimx*dimy*dimz, 0);
		else
			multiDevMemcpy<T>(arr_cpu, arr, dimx*dimy*dimz, sliceOffset);
		
		for (int i = 0; i < dimx; i++)
		{
			for (int j = 0; j < dimy; j++)
			{
				for (int k = 0; k < dimz; k++)
					printf("%12f  ", arr_cpu[i*dimy*dimz + j*dimz + k]);
				printf("\n");
			}
			printf("\n");
		}
		/**/
		delete [] arr_cpu;
	}

	template <class T>
	void printEllements(T* arr, int dimy, int dimz, int sliceOffset = 0, bool printWithOffset = false)
	{
		GPUplan* pGPUplan = GPUplan::Instance();
		int dimx = pGPUplan->getLength1D();
		dimx += (printWithOffset)?2*pGPUplan->size():0;
		
		for (int i = 0; i < dimx; i++)
		{
			for (int j = 0; j < dimy; j++)
			{
				for (int k = 0; k < dimz; k++)
					printf("%f  ", arr[i*dimy*dimz + j*dimz + k]);
				printf("\n");
			}
			printf("\n");
		}
	}

	template <typename T>
	class HostArray
	{
		T* arr;
		size_t num_elems;
	public:
		HostArray(size_t _num_elems)
		{
			num_elems = _num_elems;
			arr = new T [num_elems];
		}

		~HostArray() { delete [] arr; }

		void genRandom(){ fillRandom<T>(arr, num_elems); }

		size_t length() { return num_elems; }

		T sum() { return sumEllements(arr, num_elems); }

		void clear(T value)
		{
			for (size_t i = 0; i < num_elems; i++)
				arr[i] = value;
		}

		T* getArray() { return arr; }
	};

	template <typename T>
	class DeviceArray
	{
		T** dev_arr;
		size_t num_elems;
	public:
		DeviceArray() : dev_arr(NULL), num_elems(0){};

		DeviceArray(size_t _num_elems)
		{
			num_elems = _num_elems;
			GPUplan* pGPUplan = GPUplan::Instance();
			multiDevAlloc(dev_arr, num_elems);
			//multiDevAlloc(dev_arr, num_elems * pGPUplan->size()); // do whole hst array allocation on each node for now
		}

		~DeviceArray() { multiDevFree(dev_arr); }

		size_t length() { return num_elems; }

		T sum() { return sumEllementsMultiGPU(dev_arr, num_elems); }

		T** getArray() { return dev_arr; }
	};
}
