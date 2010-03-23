#pragma once

#include "Grid3D.h"

namespace FluidSolver3D
{
	struct ScalarField3D
	{
		int dimx, dimy, dimz;
		double dx, dy, dz;

		// access element
		double& elem(int i, int j, int k)
		{
			return u[i * dimy * dimz + j * dimz + k];
		}

		// get derivatives
		inline double d_x(int i, int j, int k)	{ return (elem(i+1, j, k) - elem(i-1, j, k)) / (2 * dx); }
		inline double d_y(int i, int j, int k)	{ return (elem(i, j+1, k) - elem(i, j-1, k)) / (2 * dy); }
		inline double d_z(int i, int j, int k)	{ return (elem(i, j, k+1) - elem(i, j, k-1)) / (2 * dz); }
		inline double d_xx(int i, int j, int k)	{ return (elem(i+1, j, k) - 2 * elem(i, j, k) + elem(i-1, j, k)) / (dx * dx); }
		inline double d_yy(int i, int j, int k)	{ return (elem(i, j+1, k) - 2 * elem(i, j, k) + elem(i, j-1, k)) / (dy * dy); }
		inline double d_zz(int i, int j, int k)	{ return (elem(i, j, k+1) - 2 * elem(i, j, k) + elem(i, j, k-1)) / (dz * dz); }

		ScalarField3D(int _dimx, int _dimy, int _dimz, double _dx, double _dy, double _dz) : 
			dimx(_dimx), dimy(_dimy), dimz(_dimz),
			dx(_dx), dy(_dy), dz(_dz)
		{
			u = new double[dimx * dimy * dimz];
		}

		ScalarField3D()
		{
			delete [] u;
		}

		void CopyFieldTo(Grid3D *grid, ScalarField3D *dest, NodeType type)
		{
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					for (int k = 0; k < dimz-1; k++)
						if (grid->GetType(i, j, k) == type)
							dest->elem(i, j, k) = elem(i, j, k);
		}

		void MergeFieldTo(Grid3D *grid, ScalarField3D *dest, NodeType type)
		{
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					for (int k = 0; k < dimz-1; k++)
						if (grid->GetType(i, j, k) == type)
							dest->elem(i, j, k) = (dest->elem(i, j, k) + elem(i, j, k)) / 2;
		}

	private:
		double *u;
	};

	struct TimeLayer3D
	{
		int dimx, dimy, dimz;
		double dx, dy, dz;

		ScalarField3D *U, *V, *W, *T;
		
		double DissFuncX(int i, int j, int k)
		{
			// TODO: eval dissipative function
			return 0.0;
		}

		double DissFuncY(int i, int j, int k)
		{
			// TODO: eval dissipative function
			return 0.0;
		}

		double DissFuncZ(int i, int j, int k)
		{
			// TODO: eval dissipative function
			return 0.0;
		}

		double DissFunc(int i, int j, int k)
		{
			return DissFuncX(i, j, k) + DissFuncY(i, j, k) + DissFunc(i, j, k);
		}

		double EvalDivError(Grid3D *grid)
		{
			double err = 0.0;
			int count = 0;
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					for (int k = 0; k < dimz-1; k++)
						if (grid->GetType(i, j, k) == IN)
						{
							// TODO: eval error
							count++;
						}
			return err / count;
		}

		void MergeLayerTo(Grid3D *grid, TimeLayer3D *dest, NodeType type)
		{
			U->MergeFieldTo(grid, dest->U, type);
			V->MergeFieldTo(grid, dest->V, type);
			W->MergeFieldTo(grid, dest->W, type);
			T->MergeFieldTo(grid, dest->T, type);
		}

		void CopyLayerTo(Grid3D *grid, TimeLayer3D *dest)
		{
			CopyLayerTo(grid, dest, IN);
			CopyLayerTo(grid, dest, OUT);
			CopyLayerTo(grid, dest, BOUND);
		}

		void CopyLayerTo(Grid3D *grid, TimeLayer3D *dest, NodeType type)
		{			
			U->CopyFieldTo(grid, dest->U, type);
			V->CopyFieldTo(grid, dest->V, type);
			W->CopyFieldTo(grid, dest->W, type);
			T->CopyFieldTo(grid, dest->T, type);
		}

		TimeLayer3D(int _dimx, int _dimy, int _dimz, double _dx, double _dy, double _dz) : 
			dimx(_dimx), dimy(_dimy), dimz(_dimz),
			dx(_dx), dy(_dy), dz(_dz)
		{
			U = new ScalarField3D(dimx, dimy, dimz, dx, dy, dz);
			V = new ScalarField3D(dimx, dimy, dimz, dx, dy, dz);
			W = new ScalarField3D(dimx, dimy, dimz, dx, dy, dz);
			T = new ScalarField3D(dimx, dimy, dimz, dx, dy, dz);
		}
		
		~TimeLayer3D()
		{
			delete U;
			delete V;
			delete W;
			delete T;
		}
	};
}