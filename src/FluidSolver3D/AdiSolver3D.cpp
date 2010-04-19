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
		temp_local = NULL;
		next_local = NULL;
		next = NULL;

		a = NULL;
		b = NULL;
		c = NULL;
		d = NULL;
		x = NULL;
	}

	void AdiSolver3D::FreeMemory()
	{
		if (cur != NULL) delete cur;
		if (temp != NULL) delete temp;
		if (half1 != NULL) delete half1;
		if (half2 != NULL) delete half2;
		if (temp_local != NULL) delete temp_local;
		if (next_local != NULL) delete next_local;
		if (next != NULL) delete next;

		if (a != NULL) delete [] a;
		if (b != NULL) delete [] b;
		if (c != NULL) delete [] c;
		if (d != NULL) delete [] d;
		if (x != NULL) delete [] x;
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

		int n = max(dimx, max(dimy, dimz));
		a = new double[n];
		b = new double[n];
		c = new double[n];
		d = new double[n];
		x = new double[n];

		cur = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);
		half1 = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);
		half2 = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);
		next = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);

		temp = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);
		temp_local = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);
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
		temp->CopyLayerTo(grid, temp_local);

		DirType dir = list[0].dir;
		for (int it = 0; it < num_local; it++)
		{
			for (size_t s = 0; s < list.size(); s++)
			{		
				SolveSegment(dt, list[s], type_U, dir, cur, temp, temp_local, next_local);
				SolveSegment(dt, list[s], type_V, dir, cur, temp, temp_local, next_local);
				SolveSegment(dt, list[s], type_W, dir, cur, temp, temp_local, next_local);
				SolveSegment(dt, list[s], type_T, dir, cur, temp, temp_local, next_local);			
			}

			// update non-linear
			if (it == 0) next_local->CopyLayerTo(grid, temp_local, IN);
				else next_local->MergeLayerTo(grid, temp_local, IN);
		}

		temp_local->CopyLayerTo(grid, temp, IN);
		next_local->CopyLayerTo(grid, next, IN);
	}

	void AdiSolver3D::SolveSegment(double dt, Segment3D seg, VarType var, DirType dir, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *temp_local, TimeLayer3D *next_local)
	{
		int n = seg.size;
		ApplyBC0(seg.posx, seg.posy, seg.posz, var, b[0], c[0], d[0]);
		BuildMatrix(dt, seg.posx, seg.posy, seg.posz, var, dir, a, b, c, d, n, cur, temp, temp_local);
		ApplyBC1(seg.endx, seg.endy, seg.endz, var, a[n-1], b[n-1], d[n-1]);
		SolveTridiagonal(a, b, c, d, x, n);
		UpdateSegment(x, seg, var, next_local);
	}

	void AdiSolver3D::UpdateSegment(double *x, Segment3D seg, VarType var, TimeLayer3D *layer)
	{
		int i = seg.posx;
		int j = seg.posy;
		int k = seg.posz;

		for (int t = 0; t < seg.size; t++)
		{
			switch (var)
			{
			case type_U: layer->U->elem(i, j, k) = x[t]; break;
			case type_V: layer->V->elem(i, j, k) = x[t]; break;
			case type_W: layer->W->elem(i, j, k) = x[t]; break;
			case type_T: layer->T->elem(i, j, k) = x[t]; break;
			}

			switch (seg.dir)
			{
			case X: i++; break;
			case Y: j++; break;
			case Z: k++; break;
			}
		}
	}

	void AdiSolver3D::BuildMatrix(double dt, int i, int j, int k, VarType var, DirType dir, double *a, double *b, double *c, double *d, int n, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *temp_local)
	{
		double vis_dx2, vis_dy2, vis_dz2;

		switch (var)
		{
		case type_U:
		case type_V:
		case type_W:
			vis_dx2 = params.v_vis / (grid->dx * grid->dx);
			vis_dy2 = params.v_vis / (grid->dy * grid->dy);
			vis_dz2 = params.v_vis / (grid->dz * grid->dz);
			break;
		case type_T:
			vis_dx2 = params.t_vis / (grid->dx * grid->dx);
			vis_dy2 = params.t_vis / (grid->dy * grid->dy);
			vis_dz2 = params.t_vis / (grid->dz * grid->dz);
			break;
		}
		
		for (int p = 1; p < n-1; p++)
		{
			switch (dir)
			{
			case X:
				a[p] = - temp_local->U->elem(i+p, j, k) / (2 * grid->dx) - vis_dx2; 
				b[p] = 1 / dt  +  2 * vis_dx2; 
				c[p] = temp_local->U->elem(i+p, j, k) / (2 * grid->dx) - vis_dx2; 
				
				switch (var)	
				{
				case type_U: d[p] = cur->U->elem(i+p, j, k) / dt - params.v_T * temp_local->T->d_x(i+p, j, k); break;
				case type_V: d[p] = cur->V->elem(i+p, j, k) / dt; break;
				case type_W: d[p] = cur->W->elem(i+p, j, k) / dt; break;
				case type_T: d[p] = cur->T->elem(i+p, j, k) / dt + params.t_phi * temp_local->DissFuncX(i+p, j, k); break;
				}	
				break;

			case Y:
				a[p] = - temp_local->V->elem(i, j+p, k) / (2 * grid->dy) - vis_dy2; 
				b[p] = 1 / dt  +  2 * vis_dy2; 
				c[p] = temp_local->V->elem(i, j+p, k) / (2 * grid->dy) - vis_dy2; 
				
				switch (var)	
				{
				case type_U: d[p] = cur->U->elem(i, j+p, k) / dt; break;
				case type_V: d[p] = cur->V->elem(i, j+p, k) / dt - params.v_T * temp_local->T->d_y(i, j+p, k); break;
				case type_W: d[p] = cur->W->elem(i, j+p, k) / dt; break;
				case type_T: d[p] = cur->T->elem(i, j+p, k) / dt + params.t_phi * temp_local->DissFuncY(i, j+p, k); break;
				}
				break;

			case Z:
				a[p] = - temp_local->W->elem(i, j, k+p) / (2 * grid->dz) - vis_dz2; 
				b[p] = 1 / dt  +  2 * vis_dz2; 
				c[p] = temp_local->W->elem(i, j, k+p) / (2 * grid->dz) - vis_dz2; 
				
				switch (var)	
				{
				case type_U: d[p] = cur->U->elem(i, j, k+p) / dt; break;
				case type_V: d[p] = cur->V->elem(i, j, k+p) / dt; break;
				case type_W: d[p] = cur->W->elem(i, j, k+p) / dt - params.v_T * temp_local->T->d_z(i, j, k+p); break;
				case type_T: d[p] = cur->T->elem(i, j, k+p) / dt + params.t_phi * temp_local->DissFuncZ(i, j, k+p); break;
				}
				break;

			}
		}
	}

	void AdiSolver3D::ApplyBC0(int i, int j, int k, VarType var, double &b0, double &c0, double &d0)
	{
		switch (grid->GetBC(i, j, k))
		{
		case NOSLIP: 
			b0 = 1.0; 
			c0 = 0.0; 
			switch (var)
			{
			case type_U: d0 = grid->GetVel(i, j, k).x; break;
			case type_V: d0 = grid->GetVel(i, j, k).y; break;
			case type_W: d0 = grid->GetVel(i, j, k).z; break;
			case type_T: d0 = grid->GetT(i, j, k); break;
			}
			break;
		case FREE: 
			b0 = 1.0; 
			c0 = -1.0; 
			d0 = 0.0; 
			break;
		}
	}

	void AdiSolver3D::ApplyBC1(int i, int j, int k, VarType var, double &a1, double &b1, double &d1)
	{
		switch (grid->GetBC(i, j, k))
		{
		case NOSLIP: 
			a1 = 0.0; 
			b1 = 1.0; 
			switch (var)
			{
			case type_U: d1 = grid->GetVel(i, j, k).x; break;
			case type_V: d1 = grid->GetVel(i, j, k).y; break;
			case type_W: d1 = grid->GetVel(i, j, k).z; break;
			case type_T: d1 = grid->GetT(i, j, k); break;
			}
			break;
		case FREE: 
			a1 = 1.0; 
			b1 = -1.0; 
			d1 = 0.0;
			break;
		}
	}	
}