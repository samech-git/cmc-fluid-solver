#include "AdiSolver3D.h"

namespace FluidSolver3D
{
	AdiSolver3D::AdiSolver3D()
	{
		grid = NULL;

		cur = NULL;
		temp = NULL;
		half1 = NULL;
		half2 = NULL;
		next_local = NULL;
		next = NULL;
	}

	void AdiSolver3D::FreeMemory()
	{
		if (cur != NULL) delete cur;
		if (temp != NULL) delete temp;
		if (half1 != NULL) delete half1;
		if (half2 != NULL) delete half2;
		if (next_local != NULL) delete next_local;
		if (next != NULL) delete next;
	}

	AdiSolver3D::~AdiSolver3D()
	{
		FreeMemory();
	}

	void AdiSolver3D::Init(Grid3D* _grid, FluidParams &_params)
	{
		grid = _grid;

		dimx = grid->dimx;
		dimy = grid->dimy;
		dimz = grid->dimz;

		params = _params;

		cur = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);
		half1 = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);
		half2 = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);
		next = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);

		temp = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);
		next_local = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);

		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				for (int k = 0; k < dimz; k++)
					switch(grid->GetType(i, j, k))
					{
					case IN:
					case BOUND:
					case OUT:
						Vec3D velocity = grid->GetVel(i, j, k);
						cur->U->elem(i, j, k) = velocity.x;
						cur->V->elem(i, j, k) = velocity.y;
						cur->W->elem(i, j, k) = velocity.z;
						cur->T->elem(i, j, k) = grid->GetT(i, j, k);
						break;
					}
	}

	void AdiSolver3D::TimeStep(double dt, int num_global, int num_local)
	{
		CreateSegments();

		cur->CopyLayerTo(grid, next);
		cur->CopyLayerTo(grid, half1);
		cur->CopyLayerTo(grid, half2);
		cur->CopyLayerTo(grid, temp);

		// do global iterations
		int it;
		double err = next->EvalDivError(grid);

		for (it = 0; (it < num_global) || (err > ERR_THRESHOLD); it++)
		{
			// alternating directions
			SolveDirection(dt, num_local, listZ, cur, temp, half1);
			SolveDirection(dt, num_local, listY, half1, temp, half2);
			SolveDirection(dt, num_local, listX, half2, temp, next);

			err = next->EvalDivError(grid);
			
			// update non-linear parameters
			if (it == 0) next->CopyLayerTo(grid, temp, IN);
				else next->MergeLayerTo(grid, temp, IN);
			
			if (it > MAX_GLOBAL_ITERS) 
			{
				printf("\nExceeded max number of iterations (%i)\n", MAX_GLOBAL_ITERS); 
				exit(1); 
			}

			if (err > ERR_THRESHOLD)
			{
				printf("\nError is too big!\n", err);
				exit(1);
			}
		}

		ClearOutterCells();

		// output error
		printf("\rerr = %.4f,", err);

		// copy result to current layer
		next->CopyLayerTo(grid, cur);
	}

	void AdiSolver3D::CreateSegments()
	{
		for (int dir = X; dir <= Z; dir++)
		{
			vector<Segment3D> *list;
			int dim1, dim2, dim3;
			switch (dir)
			{
			case X: list = &listX; dim1 = dimx; dim2 = dimy; dim3 = dimz; break;
			case Y: list = &listY; dim1 = dimy; dim2 = dimx; dim3 = dimz; break;
			case Z: list = &listZ; dim1 = dimz; dim2 = dimx; dim3 = dimy; break;
			}
			list->clear();
			
			for (int i = 0; i < dim2; i++)
				for (int j = 0; j < dim3; j++)
				{
					Segment3D seg, new_seg;
					int state = 0, incx, incy, incz;
					switch (dir)
					{
					case X:	
						seg.posx = 0; seg.posy = i; seg.posz = j; 
						incx = 1; incy = 0; incz = 0;
						break;
					case Y: 
						seg.posx = i; seg.posy = 0; seg.posz = j; 
						incx = 0; incy = 1; incz = 0;
						break;
					case Z: 
						seg.posx = i; seg.posy = j; seg.posz = 0; 
						incx = 0; incy = 0; incz = 1; 
						break;
					}
					seg.dir = (DirType)dir;

					while ((seg.posx + incx < dimx) && (seg.posy + incy < dimy) && (seg.posz + incz < dimz))
					{
						if (grid->GetType(seg.posx + incx, seg.posy + incy, seg.posz + incz) == IN)
						{
							if (state == 0) 
								new_seg = seg;
							state = 1;
						}
						else
						{
							if (state == 1)
							{
								new_seg.endx = seg.posx + incx;
								new_seg.endy = seg.posy + incy;
								new_seg.endz = seg.posz + incz;

								new_seg.size = (new_seg.endx - new_seg.posx) + (new_seg.endy - new_seg.posy) + (new_seg.endz - new_seg.posz) + 1;

								list->push_back(new_seg);
								state = 0;
							}
						}
						
						seg.posx += incx;
						seg.posy += incy;
						seg.posz += incz;
					}
				}
		}
	}

	void AdiSolver3D::SolveDirection(double dt, int num_local, vector<Segment3D> &list, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *next)
	{
	}

	void AdiSolver3D::SolveSegment(double dt, Segment3D seg, VarType var, DirType dir, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *temp_local, TimeLayer3D *next_local)
	{
	}

	void AdiSolver3D::UpdateSegment(double *x, Segment3D seg, VarType var, TimeLayer3D *layer)
	{
	}

	void AdiSolver3D::BuildMatrix(double dt, int i, int j, int k, VarType var, DirType dir, double *a, double *b, double *c, double *d, int n, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *temp_local)
	{
	}

	void AdiSolver3D::ApplyBC0(int i, int j, int k, VarType var, double &b0, double &c0, double &d0)
	{
	}

	void AdiSolver3D::ApplyBC1(int i, int j, int k, VarType var, double &a1, double &b1, double &d1)
	{
	}	
}