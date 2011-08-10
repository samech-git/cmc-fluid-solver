#include "test_util.h"

namespace TestUtil
{

	double sumEllements(FluidSolver3D::Node* arr, size_t num_elems)
	{
		double sum = 0.;
		for (size_t i = 0; i < num_elems; i++)
			sum += arr[i].T;
		return sum;
	}

	double sumEllementsMultiGPU(FluidSolver3D::Node** arr, size_t num_elems_total)
	{
		double sum = 0.;
		FluidSolver3D::Node *arr_cpu = new FluidSolver3D::Node[num_elems_total];
		multiDevMemcpy<FluidSolver3D::Node>(arr_cpu, arr, num_elems_total);
		sum = sumEllements(arr_cpu,  num_elems_total);
		delete [] arr_cpu;
		return sum;
	}
}
