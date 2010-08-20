#pragma once

#include "Grid3D.h"

using namespace FluidSolver3D;

#include <cuda_runtime.h>

extern void CopyFieldTo_GPU(int dimx, int dimy, int dimz, FTYPE *src, FTYPE *dest, Node *nodes, NodeType target);
extern void MergeFieldTo_GPU(int dimx, int dimy, int dimz, FTYPE *src, FTYPE *dest, Node *nodes, NodeType target);
extern void CopyFromGrid_GPU(int dimx, int dimy, int dimz, FTYPE *u, FTYPE *v, FTYPE *w, FTYPE *T, Node *nodes, NodeType target);
extern void Clear_GPU(int dimx, int dimy, int dimz, FTYPE *u, FTYPE *v, FTYPE *w, FTYPE *T, Node *nodes, NodeType target, FTYPE const_u, FTYPE const_v, FTYPE const_w, FTYPE const_T);
extern void Transpose_GPU(int dimx, int dimy, int dimz, FTYPE *u, FTYPE *dest_u);

enum BackendType { CPU, GPU };

namespace FluidSolver3D
{
	struct ScalarField3D
	{
		BackendType hw; 
		int dimx, dimy, dimz;
		FTYPE dx, dy, dz;

		// access element
		inline FTYPE& elem(int i, int j, int k)
		{
			return u[i * dimy * dimz + j * dimz + k];
		}

		// access the whole array
		FTYPE *getArray()
		{
			return u;
		}

		// get derivatives
		inline FTYPE d_x(int i, int j, int k)	{ return (elem(i+1, j, k) - elem(i-1, j, k)) / (2 * dx); }
		inline FTYPE d_y(int i, int j, int k)	{ return (elem(i, j+1, k) - elem(i, j-1, k)) / (2 * dy); }
		inline FTYPE d_z(int i, int j, int k)	{ return (elem(i, j, k+1) - elem(i, j, k-1)) / (2 * dz); }
		inline FTYPE d_xx(int i, int j, int k)	{ return (elem(i+1, j, k) - 2 * elem(i, j, k) + elem(i-1, j, k)) / (dx * dx); }
		inline FTYPE d_yy(int i, int j, int k)	{ return (elem(i, j+1, k) - 2 * elem(i, j, k) + elem(i, j-1, k)) / (dy * dy); }
		inline FTYPE d_zz(int i, int j, int k)	{ return (elem(i, j, k+1) - 2 * elem(i, j, k) + elem(i, j, k-1)) / (dz * dz); }

		ScalarField3D(BackendType _hw, int _dimx, int _dimy, int _dimz, FTYPE _dx, FTYPE _dy, FTYPE _dz) : 
			hw(_hw), dimx(_dimx), dimy(_dimy), dimz(_dimz),
			dx(_dx), dy(_dy), dz(_dz)
		{
			switch( hw )
			{
			case CPU: u = new FTYPE[dimx * dimy * dimz]; break;
			case GPU: cudaMalloc(&u, dimx * dimy * dimz * sizeof(FTYPE)); break;
			}
		}

		ScalarField3D(BackendType _hw, ScalarField3D *field) : 
			hw(_hw), dimx(field->dimx), dimy(field->dimy), dimz(field->dimz),
			dx(field->dx), dy(field->dy), dz(field->dz)
		{
			switch( hw )
			{
			case CPU: 
				u = new FTYPE[dimx * dimy * dimz]; 
				switch( field->hw )
				{
				case CPU: cudaMemcpy(u, field->getArray(), dimx * dimy * dimz * sizeof(FTYPE), cudaMemcpyHostToHost); break;
				case GPU: cudaMemcpy(u, field->getArray(), dimx * dimy * dimz * sizeof(FTYPE), cudaMemcpyDeviceToHost); break;
				}
				break;
			case GPU: 
				cudaMalloc(&u, dimx * dimy * dimz * sizeof(FTYPE)); 
				switch( field->hw )
				{
				case CPU: cudaMemcpy(u, field->getArray(), dimx * dimy * dimz * sizeof(FTYPE), cudaMemcpyHostToDevice); break;
				case GPU: cudaMemcpy(u, field->getArray(), dimx * dimy * dimz * sizeof(FTYPE), cudaMemcpyDeviceToDevice); break;
				}
				break;
			}
		}

		~ScalarField3D()
		{
			switch( hw )
			{
			case CPU: delete [] u; break;
			case GPU: cudaFree(u); break;
			}
		}

		void CopyFieldTo(Grid3D *grid, ScalarField3D *dest, NodeType type)
		{
			switch( hw )
			{
			case CPU:
				{
					for (int i = 0; i < dimx; i++)
						for (int j = 0; j < dimy; j++)
							for (int k = 0; k < dimz; k++)
								if (grid->GetType(i, j, k) == type)
									dest->elem(i, j, k) = elem(i, j, k);
					break;
				}
			case GPU: 
				{
					CopyFieldTo_GPU(dimx, dimy, dimz, u, dest->getArray(), grid->GetNodesGPU(), type); 
					break;
				}
			}
		}

		void MergeFieldTo(Grid3D *grid, ScalarField3D *dest, NodeType type)
		{
			switch( hw )
			{
			case CPU:
				{
					#pragma omp parallel default(none) firstprivate(type) shared(grid, dest)
					{
						#pragma omp for
						for (int i = 0; i < dimx; i++)
							for (int j = 0; j < dimy; j++)
								for (int k = 0; k < dimz; k++)
									if (grid->GetType(i, j, k) == type)
										dest->elem(i, j, k) = (dest->elem(i, j, k) + elem(i, j, k)) / 2;
					}
					break;
				}
			case GPU:
				{
					MergeFieldTo_GPU(dimx, dimy, dimz, u, dest->getArray(), grid->GetNodesGPU(), type);
					break;
				}
			}
		}

		void Transpose(ScalarField3D *dest)
		{
			switch( hw )
			{
			case CPU:
				{
					printf("Transpose is not implemented for CPU!\n");
				}
			case GPU:
				{
					Transpose_GPU(dimx, dimy, dimz, getArray(), dest->getArray());
					break;	
				}
			}		
		}

	private:
		FTYPE *u;		// store on CPU or GPU according to the hw flag
	};

	struct TimeLayer3D
	{
		BackendType hw;
		int dimx, dimy, dimz;
		FTYPE dx, dy, dz;

		ScalarField3D *U, *V, *W, *T;
		
		inline FTYPE DissFuncX(int i, int j, int k)
		{
			FTYPE u_x = U->d_x(i, j, k);
			FTYPE v_x = V->d_x(i, j, k);
			FTYPE w_x = W->d_x(i, j, k);

			FTYPE u_y = U->d_y(i, j, k);
			FTYPE u_z = U->d_z(i, j, k);

			return 2 * u_x * u_x + v_x * v_x + w_x * w_x + v_x * u_y + w_x * u_z;
		}

		inline FTYPE DissFuncY(int i, int j, int k)
		{
			FTYPE u_y = U->d_y(i, j, k);
			FTYPE v_y = V->d_y(i, j, k);
			FTYPE w_y = W->d_y(i, j, k);

			FTYPE v_x = V->d_x(i, j, k);
			FTYPE v_z = V->d_z(i, j, k);

			return u_y * u_y + 2 * v_y * v_y + w_y * w_y + u_y * v_x + w_y * v_z;
		}

		inline FTYPE DissFuncZ(int i, int j, int k)
		{
			FTYPE u_z = U->d_z(i, j, k);
			FTYPE v_z = V->d_z(i, j, k);
			FTYPE w_z = W->d_z(i, j, k);

			FTYPE w_x = W->d_x(i, j, k);
			FTYPE w_y = W->d_y(i, j, k);

			return u_z * u_z + v_z * v_z + 2 * w_z * w_z + u_z * w_x + v_z * w_y;
		}

		inline FTYPE DissFunc(int i, int j, int k)
		{
			return DissFuncX(i, j, k) + DissFuncY(i, j, k) + DissFunc(i, j, k);
		}

		double EvalDivError(Grid3D *grid)
		{
			ScalarField3D *U_cpu = new ScalarField3D(CPU, U);
			ScalarField3D *V_cpu = new ScalarField3D(CPU, V);
			ScalarField3D *W_cpu = new ScalarField3D(CPU, W);

			double err = 0.0;
			int count = 0;
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					for (int k = 0; k < dimz-1; k++)
						if (grid->GetType(i, j, k) == NODE_IN)
						{
							double err_x = (U_cpu->elem(i, j, k) + U_cpu->elem(i, j-1, k) + U_cpu->elem(i, j-1, k-1) + U_cpu->elem(i, j, k-1) -
								U_cpu->elem(i-1, j, k) - U_cpu->elem(i-1, j-1, k) - U_cpu->elem(i-1, j-1, k-1) - U_cpu->elem(i-1, j, k-1)) * dz * dy / 4.0;

							double err_y = (V_cpu->elem(i, j, k) + V_cpu->elem(i-1, j, k) + V_cpu->elem(i-1, j, k-1) + V_cpu->elem(i, j, k-1) -
								V_cpu->elem(i, j-1, k) - V_cpu->elem(i-1, j-1, k) - V_cpu->elem(i-1, j-1, k-1) - V_cpu->elem(i, j-1, k-1)) * dx * dz / 4.0;

							double err_z = (W_cpu->elem(i, j, k) + W_cpu->elem(i, j-1, k) + W_cpu->elem(i-1, j-1, k) + W_cpu->elem(i-1, j, k) -
								W_cpu->elem(i, j, k-1) - W_cpu->elem(i, j-1, k-1) - W_cpu->elem(i-1, j-1, k-1) - W_cpu->elem(i-1, j, k-1)) * dx * dy / 4.0;

							err += abs(err_x + err_y + err_z);
							count++;
						}

			delete U_cpu;
			delete V_cpu;
			delete W_cpu;

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
			cudaMemcpyKind dir;
			switch( hw )
			{
			case CPU:
				switch( dest->hw )
				{
				case CPU: dir = cudaMemcpyHostToHost; break;
				case GPU: dir = cudaMemcpyHostToDevice; break;
				}
				break;
			case GPU:
				switch( dest->hw )
				{
				case CPU: dir = cudaMemcpyDeviceToHost; break;
				case GPU: dir = cudaMemcpyDeviceToDevice; break;
				}
				break;
			}

			cudaMemcpy(dest->U->getArray(), U->getArray(), dimx * dimy * dimz * sizeof(FTYPE), dir);
			cudaMemcpy(dest->V->getArray(), V->getArray(), dimx * dimy * dimz * sizeof(FTYPE), dir);
			cudaMemcpy(dest->W->getArray(), W->getArray(), dimx * dimy * dimz * sizeof(FTYPE), dir);
			cudaMemcpy(dest->T->getArray(), T->getArray(), dimx * dimy * dimz * sizeof(FTYPE), dir);
		}

		void CopyLayerTo(Grid3D *grid, TimeLayer3D *dest, NodeType type)
		{			
			U->CopyFieldTo(grid, dest->U, type);
			V->CopyFieldTo(grid, dest->V, type);
			W->CopyFieldTo(grid, dest->W, type);
			T->CopyFieldTo(grid, dest->T, type);
		}

		void CopyFromGrid(Grid3D *grid)
		{
			switch( hw )
			{
			case CPU: 
				{
					for (int i = 0; i < dimx; i++)
						for (int j = 0; j < dimy; j++)
							for (int k = 0; k < dimz; k++)
							{
								Vec3D vel = grid->GetVel(i, j, k);
								U->elem(i, j, k) = vel.x;
								V->elem(i, j, k) = vel.y;
								W->elem(i, j, k) = vel.z;
								T->elem(i, j, k) = (FTYPE)grid->GetT(i, j, k);
							}
					break;
				}
			case GPU: 
				{
					FTYPE *u_cpu = new FTYPE[dimx * dimy * dimz];
					FTYPE *v_cpu = new FTYPE[dimx * dimy * dimz];
					FTYPE *w_cpu = new FTYPE[dimx * dimy * dimz];
					FTYPE *T_cpu = new FTYPE[dimx * dimy * dimz];

					for (int i = 0; i < dimx; i++)
						for (int j = 0; j < dimy; j++)
							for (int k = 0; k < dimz; k++)
							{
								int id = i * dimy * dimz + j * dimz + k;
								Vec3D vel = grid->GetVel(i, j, k);
								u_cpu[id] = vel.x;
								v_cpu[id] = vel.y;
								w_cpu[id] = vel.z;
								T_cpu[id] = grid->GetT(i, j, k);
							}


					cudaMemcpy(U->getArray(), u_cpu, sizeof(FTYPE) * dimx * dimy * dimz, cudaMemcpyHostToDevice);
					cudaMemcpy(V->getArray(), v_cpu, sizeof(FTYPE) * dimx * dimy * dimz, cudaMemcpyHostToDevice);
					cudaMemcpy(W->getArray(), w_cpu, sizeof(FTYPE) * dimx * dimy * dimz, cudaMemcpyHostToDevice);
					cudaMemcpy(T->getArray(), T_cpu, sizeof(FTYPE) * dimx * dimy * dimz, cudaMemcpyHostToDevice);

					delete [] u_cpu;
					delete [] v_cpu;
					delete [] w_cpu;
					delete [] T_cpu;
					break;
				}
			}
		}

		void CopyToGrid(Grid3D *grid)
		{
			switch( hw )
			{
			case CPU:
				{
					for (int i = 0; i < dimx; i++)
						for (int j = 0; j < dimy; j++)
							for (int k = 0; k < dimz; k++)
								grid->SetNodeVel(i, j, k, Vec3D(U->elem(i, j, k), V->elem(i, j, k), W->elem(i, j , k)));
					break;
				}
			case GPU:
				{
					FTYPE *u_cpu = new FTYPE[dimx * dimy * dimz];
					FTYPE *v_cpu = new FTYPE[dimx * dimy * dimz];
					FTYPE *w_cpu = new FTYPE[dimx * dimy * dimz];
					cudaMemcpy(u_cpu, U->getArray(), sizeof(FTYPE) * dimx * dimy * dimz, cudaMemcpyDeviceToHost);
					cudaMemcpy(v_cpu, V->getArray(), sizeof(FTYPE) * dimx * dimy * dimz, cudaMemcpyDeviceToHost);
					cudaMemcpy(w_cpu, W->getArray(), sizeof(FTYPE) * dimx * dimy * dimz, cudaMemcpyDeviceToHost);

					for (int i = 0; i < dimx; i++)
						for (int j = 0; j < dimy; j++)
							for (int k = 0; k < dimz; k++)
								grid->SetNodeVel(i, j, k, Vec3D(u_cpu[i * dimy * dimz + j * dimz + k],
																v_cpu[i * dimy * dimz + j * dimz + k],
																w_cpu[i * dimy * dimz + j * dimz + k]));

					delete [] u_cpu;
					delete [] v_cpu;
					delete [] w_cpu;

					break;
				}
			}
		}

		void FilterToArrays(Vec3D *outV, double *outT, int outdimx, int outdimy, int outdimz)
		{
			if (outdimx == 0) outdimx = dimx;
			if (outdimy == 0) outdimy = dimy;
			if (outdimz == 0) outdimz = dimz;

			switch( hw )
			{
			case CPU:
				{
					for (int i = 0; i < outdimx; i++)
						for (int j = 0; j < outdimy; j++)
							for (int k = 0; k < outdimz; k++)
							{	
								int x = (i * dimx / outdimx);
								int y = (j * dimy / outdimy);
								int z = (k * dimz / outdimz);
								int ind = i * outdimy * outdimz + j * outdimz + k;
								outV[ind].x = U->elem(x, y, z); 
								outV[ind].y = V->elem(x, y, z);
								outV[ind].z = W->elem(x, y, z);
								outT[ind] = T->elem(x, y, z);
							}
					break;
				}
			case GPU:
				{
					FTYPE *u_cpu = new FTYPE[dimx * dimy * dimz];
					FTYPE *v_cpu = new FTYPE[dimx * dimy * dimz];
					FTYPE *w_cpu = new FTYPE[dimx * dimy * dimz];
					FTYPE *T_cpu = new FTYPE[dimx * dimy * dimz];
					cudaMemcpy(u_cpu, U->getArray(), sizeof(FTYPE) * dimx * dimy * dimz, cudaMemcpyDeviceToHost);
					cudaMemcpy(v_cpu, V->getArray(), sizeof(FTYPE) * dimx * dimy * dimz, cudaMemcpyDeviceToHost);
					cudaMemcpy(w_cpu, W->getArray(), sizeof(FTYPE) * dimx * dimy * dimz, cudaMemcpyDeviceToHost);
					cudaMemcpy(T_cpu, T->getArray(), sizeof(FTYPE) * dimx * dimy * dimz, cudaMemcpyDeviceToHost);

					for (int i = 0; i < outdimx; i++)
						for (int j = 0; j < outdimy; j++)
							for (int k = 0; k < outdimz; k++)
							{	
								int x = (i * dimx / outdimx);
								int y = (j * dimy / outdimy);
								int z = (k * dimz / outdimz);
								int ind = i * outdimy * outdimz + j * outdimz + k;
								outV[ind].x = u_cpu[x * dimy * dimz + y * dimz + z]; 
								outV[ind].y = v_cpu[x * dimy * dimz + y * dimz + z]; 
								outV[ind].z = w_cpu[x * dimy * dimz + y * dimz + z]; 
								outT[ind] = T_cpu[x * dimy * dimz + y * dimz + z]; 
							}

					delete [] u_cpu;
					delete [] v_cpu;
					delete [] w_cpu;
					delete [] T_cpu;
					break;
				}
			}
		}

		void CopyFromGrid(Grid3D *grid, NodeType target)
		{
			switch( hw )
			{
			case CPU:
				{
					for (int i = 0; i < dimx; i++)
						for (int j = 0; j < dimy; j++)
							for (int k = 0; k < dimz; k++)
								if (grid->GetType(i, j, k) == target)
								{
									Vec3D velocity = grid->GetVel(i, j, k);
									U->elem(i, j, k) = (FTYPE)velocity.x;
									V->elem(i, j, k) = (FTYPE)velocity.y;
									W->elem(i, j, k) = (FTYPE)velocity.z;
									T->elem(i, j, k) = (FTYPE)grid->GetT(i, j, k);
								}
					break;
				}
			case GPU:
				{
					CopyFromGrid_GPU(dimx, dimy, dimz, U->getArray(), V->getArray(), W->getArray(), T->getArray(), grid->GetNodesGPU(), target); 
					break;
				}
			}
		}

		void Clear(Grid3D *grid, NodeType target, FTYPE const_u, FTYPE const_v, FTYPE const_w, FTYPE const_T)
		{
			switch( hw )
			{
			case CPU:
				{
					for (int i = 0; i < dimx; i++)
						for (int j = 0; j < dimy; j++)
							for (int k = 0; k < dimz; k++)
								if (grid->GetType(i, j, k) == target)
								{
									U->elem(i, j, k) = const_u;
									V->elem(i, j, k) = const_v;
									W->elem(i, j, k) = const_w;
									T->elem(i, j, k) = const_T;
								}
					break;
				}
			case GPU:
				{
					Clear_GPU(dimx, dimy, dimz, U->getArray(), V->getArray(), W->getArray(), T->getArray(), grid->GetNodesGPU(), target, const_u, const_v, const_w, const_T);
					break;
				}
			}
		}

		void Transpose(TimeLayer3D *dest)
		{
			U->Transpose(dest->U);
			V->Transpose(dest->V);
			W->Transpose(dest->W);
			T->Transpose(dest->T);
		}

		void PrintToFile(FILE *file, const char *desc, FTYPE *u, int dimx, int dimy, int dimz)
		{
			fprintf(file, "Array %s[%i,%i,%i]:\n", desc, dimx, dimy, dimz);
			for (int i = 0; i < dimx; i++)
			{			
				for (int j = 0; j < dimy; j++)
				{
					for (int k = 0; k < dimz; k++)
					{
						int id = i * dimy * dimz + j * dimz + k;
						fprintf(file, "%.8f ", u[id]);
					}
					fprintf(file, "\n");
				}
				fprintf(file, "\n");
			}
		}

		void PrintToFile(const char *filename)
		{
			FILE *file = NULL;
			fopen_s(&file, filename, "w");

			if( hw == GPU )
			{
				FTYPE *u_cpu = new FTYPE[dimx * dimy * dimz];
				FTYPE *v_cpu = new FTYPE[dimx * dimy * dimz];
				FTYPE *w_cpu = new FTYPE[dimx * dimy * dimz];
				FTYPE *T_cpu = new FTYPE[dimx * dimy * dimz];
				cudaMemcpy(u_cpu, U->getArray(), sizeof(FTYPE) * dimx * dimy * dimz, cudaMemcpyDeviceToHost);
				cudaMemcpy(v_cpu, V->getArray(), sizeof(FTYPE) * dimx * dimy * dimz, cudaMemcpyDeviceToHost);
				cudaMemcpy(w_cpu, W->getArray(), sizeof(FTYPE) * dimx * dimy * dimz, cudaMemcpyDeviceToHost);
				cudaMemcpy(T_cpu, T->getArray(), sizeof(FTYPE) * dimx * dimy * dimz, cudaMemcpyDeviceToHost);

				PrintToFile(file, "u", u_cpu, dimx, dimy, dimz);
				PrintToFile(file, "v", v_cpu, dimx, dimy, dimz);
				PrintToFile(file, "w", w_cpu, dimx, dimy, dimz);
				PrintToFile(file, "T", T_cpu, dimx, dimy, dimz);
				
				delete [] u_cpu;
				delete [] v_cpu;
				delete [] w_cpu;
				delete [] T_cpu;
			}
			else
				printf("Print to file is not implemented for CPU! :)\n");

			fclose(file);
		}
		
		TimeLayer3D(BackendType _hw, int _dimx, int _dimy, int _dimz, FTYPE _dx, FTYPE _dy, FTYPE _dz) : 
			hw(_hw), dimx(_dimx), dimy(_dimy), dimz(_dimz),
			dx(_dx), dy(_dy), dz(_dz)
		{
			U = new ScalarField3D(hw, dimx, dimy, dimz, dx, dy, dz);
			V = new ScalarField3D(hw, dimx, dimy, dimz, dx, dy, dz);
			W = new ScalarField3D(hw, dimx, dimy, dimz, dx, dy, dz);
			T = new ScalarField3D(hw, dimx, dimy, dimz, dx, dy, dz);
		}

		TimeLayer3D(BackendType _hw, Grid3D *grid) : 
			hw(_hw), dimx(grid->dimx), dimy(grid->dimy), dimz(grid->dimz),
			dx((FTYPE)grid->dx), dy((FTYPE)grid->dy), dz((FTYPE)grid->dz)
		{
			U = new ScalarField3D(hw, dimx, dimy, dimz, dx, dy, dz);
			V = new ScalarField3D(hw, dimx, dimy, dimz, dx, dy, dz);
			W = new ScalarField3D(hw, dimx, dimy, dimz, dx, dy, dz);
			T = new ScalarField3D(hw, dimx, dimy, dimz, dx, dy, dz);

			CopyFromGrid(grid);
		}
		
		~TimeLayer3D()
		{
			delete U;
			delete V;
			delete W;
			delete T;
		}
	};

	struct TimeLayer3D_GPU
	{
		FTYPE *u, *v, *w, *T;
		int dimx, dimy, dimz;

		inline __device__ FTYPE& elem(FTYPE *arr, int i, int j, int k)
		{
			return arr[i * dimy * dimz + j * dimz + k];
		}

		inline __device__ FTYPE d_x(FTYPE *arr, FTYPE dx, int i, int j, int k)	{ return (elem(arr, i+1, j, k) - elem(arr, i-1, j, k)) / (2 * dx); }
		inline __device__ FTYPE d_y(FTYPE *arr, FTYPE dy, int i, int j, int k)	{ return (elem(arr, i, j+1, k) - elem(arr, i, j-1, k)) / (2 * dy); }
		inline __device__ FTYPE d_z(FTYPE *arr, FTYPE dz, int i, int j, int k)	{ return (elem(arr, i, j, k+1) - elem(arr, i, j, k-1)) / (2 * dz); }
		inline __device__ FTYPE d_xx(FTYPE *arr, FTYPE dx, int i, int j, int k)	{ return (elem(arr, i+1, j, k) - 2 * elem(arr, i, j, k) + elem(arr, i-1, j, k)) / (dx * dx); }
		inline __device__ FTYPE d_yy(FTYPE *arr, FTYPE dy, int i, int j, int k)	{ return (elem(arr, i, j+1, k) - 2 * elem(arr, i, j, k) + elem(arr, i, j-1, k)) / (dy * dy); }
		inline __device__ FTYPE d_zz(FTYPE *arr, FTYPE dz, int i, int j, int k)	{ return (elem(arr, i, j, k+1) - 2 * elem(arr, i, j, k) + elem(arr, i, j, k-1)) / (dz * dz); }

		inline __device__ FTYPE DissFuncX(FTYPE dx, FTYPE dy, FTYPE dz, int i, int j, int k)
		{
			FTYPE u_x = d_x(u, dx, i, j, k);
			FTYPE v_x = d_x(v, dx, i, j, k);
			FTYPE w_x = d_x(w, dx, i, j, k);

			FTYPE u_y = d_y(u, dy, i, j, k);
			FTYPE u_z = d_z(u, dz, i, j, k);

			return 2 * u_x * u_x + v_x * v_x + w_x * w_x + v_x * u_y + w_x * u_z;
		}

		inline __device__ FTYPE DissFuncY(FTYPE dx, FTYPE dy, FTYPE dz, int i, int j, int k)
		{
			FTYPE u_y = d_y(u, dy, i, j, k);
			FTYPE v_y = d_y(v, dy, i, j, k);
			FTYPE w_y = d_y(w, dy, i, j, k);

			FTYPE v_x = d_x(v, dx, i, j, k);
			FTYPE v_z = d_z(v, dz, i, j, k);

			return u_y * u_y + 2 * v_y * v_y + w_y * w_y + u_y * v_x + w_y * v_z;
		}

		inline __device__ FTYPE DissFuncZ(FTYPE dx, FTYPE dy, FTYPE dz, int i, int j, int k)
		{
			FTYPE u_z = d_z(u, dz, i, j, k);
			FTYPE v_z = d_z(v, dz, i, j, k);
			FTYPE w_z = d_z(w, dz, i, j, k);

			FTYPE w_x = d_x(w, dx, i, j, k);
			FTYPE w_y = d_y(w, dy, i, j, k);

			return u_z * u_z + v_z * v_z + 2 * w_z * w_z + u_z * w_x + v_z * w_y;
		}

		inline __device__ FTYPE DissFuncZ_as_Y(FTYPE dx, FTYPE dy, FTYPE dz, int i, int j, int k)
		{
			FTYPE u_z = d_y(u, dy, i, j, k);
			FTYPE v_z = d_y(v, dy, i, j, k);
			FTYPE w_z = d_y(w, dy, i, j, k);

			FTYPE w_x = d_x(w, dx, i, j, k);
			FTYPE w_y = d_z(w, dz, i, j, k);

			return u_z * u_z + v_z * v_z + 2 * w_z * w_z + u_z * w_x + v_z * w_y;
		}

		TimeLayer3D_GPU(TimeLayer3D *layer)
		{
			u = layer->U->getArray();
			v = layer->V->getArray();
			w = layer->W->getArray();
			T = layer->T->getArray();

			dimx = layer->dimx;
			dimy = layer->dimy;
			dimz = layer->dimz;
		}
	};
}