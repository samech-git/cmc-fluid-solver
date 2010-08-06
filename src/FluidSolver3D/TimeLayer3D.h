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
			#pragma omp parallel default(none) firstprivate(type) shared(grid, dest)
			{
				#pragma omp for
				for (int i = 0; i < dimx-1; i++)
					for (int j = 0; j < dimy-1; j++)
						for (int k = 0; k < dimz-1; k++)
							if (grid->GetType(i, j, k) == type)
								dest->elem(i, j, k) = (dest->elem(i, j, k) + elem(i, j, k)) / 2;
			}
		}

		void Print(char *filename)
		{
			FILE *file = NULL;
			fopen_s(&file, filename, "w");
			for (int k = 0; k < dimz-1; k++)
			{
				for (int i = 0; i < dimx-1; i++)
				{
					for (int j = 0; j < dimy-1; j++)
						fprintf(file, "%.2f ", elem(i, j, k));
					fprintf(file, "\n");
				}
				fprintf(file, "\n");
			}
			fclose(file);
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
			double u_x = U->d_x(i, j, k);
			double v_x = V->d_x(i, j, k);
			double w_x = W->d_x(i, j, k);

			double u_y = U->d_y(i, j, k);
			double u_z = U->d_z(i, j, k);

			return 2 * u_x * u_x + v_x * v_x + w_x * w_x + v_x * u_y + w_x * u_z;
		}

		double DissFuncY(int i, int j, int k)
		{
			double u_y = U->d_y(i, j, k);
			double v_y = V->d_y(i, j, k);
			double w_y = W->d_y(i, j, k);

			double v_x = V->d_x(i, j, k);
			double v_z = V->d_z(i, j, k);

			return u_y * u_y + 2 * v_y * v_y + w_y * w_y + u_y * v_x + w_y * v_z;
		}

		double DissFuncZ(int i, int j, int k)
		{
			double u_z = U->d_z(i, j, k);
			double v_z = V->d_z(i, j, k);
			double w_z = W->d_z(i, j, k);

			double w_x = W->d_x(i, j, k);
			double w_y = W->d_y(i, j, k);

			return u_z * u_z + v_z * v_z + 2 * w_z * w_z + u_z * w_x + v_z * w_y;
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
						if (grid->GetType(i, j, k) == NODE_IN)
						{
							double err_x = (U->elem(i, j, k) + U->elem(i, j-1, k) + U->elem(i, j-1, k-1) + U->elem(i, j, k-1) -
								U->elem(i-1, j, k) - U->elem(i-1, j-1, k) - U->elem(i-1, j-1, k-1) - U->elem(i-1, j, k-1)) * dz * dy / 4.0;

							double err_y = (V->elem(i, j, k) + V->elem(i-1, j, k) + V->elem(i-1, j, k-1) + V->elem(i, j, k-1) -
								V->elem(i, j-1, k) - V->elem(i-1, j-1, k) - V->elem(i-1, j-1, k-1) - V->elem(i, j-1, k-1)) * dx * dz / 4.0;

							double err_z = (W->elem(i, j, k) + W->elem(i, j-1, k) + W->elem(i-1, j-1, k) + W->elem(i-1, j, k) -
								W->elem(i, j, k-1) - W->elem(i, j-1, k-1) - W->elem(i-1, j-1, k-1) - W->elem(i-1, j, k-1)) * dx * dy / 4.0;

							err += abs(err_x + err_y + err_z);
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

		void CopyLayerTo(TimeLayer3D *dest)
		{
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					for (int k = 0; k < dimz-1; k++)
					{
						dest->U->elem(i, j, k) = U->elem(i, j, k);
						dest->V->elem(i, j, k) = V->elem(i, j, k);
						dest->W->elem(i, j, k) = W->elem(i, j, k);
						dest->T->elem(i, j, k) = T->elem(i, j, k);
					}
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