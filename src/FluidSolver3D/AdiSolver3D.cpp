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
		if (next != NULL) delete next;

		if (a != NULL) delete [] a;
		if (b != NULL) delete [] b;
		if (c != NULL) delete [] c;
		if (d != NULL) delete [] d;
		if (x != NULL) delete [] x;
	}

	AdiSolver3D::~AdiSolver3D()
	{
		prof.PrintTimings();

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
		a = new double[n * n * n];
		b = new double[n * n * n];
		c = new double[n * n * n];
		d = new double[n * n * n];
		x = new double[n * n * n];

		cur = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);
		half1 = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);
		half2 = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);
		next = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);
		temp = new TimeLayer3D(grid->dimx, grid->dimy, grid->dimz, grid->dx, grid->dy, grid->dz);

		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				for (int k = 0; k < dimz; k++)
					switch(grid->GetType(i, j, k))
					{
					case NODE_IN:
					case NODE_BOUND:
					case NODE_VALVE:
					case NODE_OUT:
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

		prof.StartEvent();

		cur->CopyLayerTo(temp);

		prof.StopEvent("CopyLayer");

		// do global iterations
		for (int it = 0; it < num_global; it++)
		{
			// alternating directions
			SolveDirection(dt, num_local, listZ, cur, temp, half1);		
			SolveDirection(dt, num_local, listY, half1, temp, half2);
			SolveDirection(dt, num_local, listX, half2, temp, next);

			prof.StartEvent();
				
			// update non-linear parameters
			next->MergeLayerTo(grid, temp, NODE_IN);

			prof.StopEvent("MergeLayer");
		}

		prof.StartEvent();

		// output error
		double err = next->EvalDivError(grid);

		prof.StopEvent("EvalDivError");

		if (err > ERR_THRESHOLD) {
			printf("\nError is too big!\n", err);
			exit(1);
		}
		else
			printf("\rerr = %.4f,", err);

		prof.StartEvent();

		ClearOutterCells();

		prof.StopEvent("ClearLayer");

		// swap current/next pointers 
		TimeLayer3D *tmp = next;
		next = cur;
		cur = tmp;
	}

	void AdiSolver3D::CreateSegments()
	{
		prof.StartEvent();

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
						if (grid->GetType(seg.posx + incx, seg.posy + incy, seg.posz + incz) == NODE_IN)
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

		prof.StopEvent("CreateSegments");
	}

	void AdiSolver3D::SolveDirection(double dt, int num_local, vector<Segment3D> &list, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *next)
	{
		DirType dir = list[0].dir;
		for (int it = 0; it < num_local; it++)
		{
			prof.StartEvent();

			#pragma omp parallel default(none) firstprivate(dt, dir) shared(list, cur, temp, next)
			{
				#pragma omp for
				for (int s = 0; s < (int)list.size(); s++)
				{		
					SolveSegment(dt, s, list[s], type_U, dir, cur, temp, next);
					SolveSegment(dt, s, list[s], type_V, dir, cur, temp, next);
					SolveSegment(dt, s, list[s], type_W, dir, cur, temp, next);
					SolveSegment(dt, s, list[s], type_T, dir, cur, temp, next);			
				}
			}

			switch( dir )
			{
			case X: prof.StopEvent("SolveSegments_X"); break;
			case Y: prof.StopEvent("SolveSegments_Y"); break;
			case Z: prof.StopEvent("SolveSegments_Z"); break;
			}

			prof.StartEvent();

			// update non-linear
			next->MergeLayerTo(grid, temp, NODE_IN);

			prof.StopEvent("MergeLayer");
		}
	}

	void AdiSolver3D::SolveSegment(double dt, int id, Segment3D seg, VarType var, DirType dir, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *next)
	{
		int n = seg.size;

		int max_n = max(dimx, max(dimy, dimz));
		double *a = this->a + id * max_n;
		double *b = this->b + id * max_n;
		double *c = this->c + id * max_n;
		double *d = this->d + id * max_n;
		double *x = this->x + id * max_n;

		ApplyBC0(seg.posx, seg.posy, seg.posz, var, b[0], c[0], d[0]);
		ApplyBC1(seg.endx, seg.endy, seg.endz, var, a[n-1], b[n-1], d[n-1]);
		BuildMatrix(dt, seg.posx, seg.posy, seg.posz, var, dir, a, b, c, d, n, cur, temp);
		
		SolveTridiagonal(a, b, c, d, x, n);
		
		UpdateSegment(x, seg, var, next);
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

	void AdiSolver3D::BuildMatrix(double dt, int i, int j, int k, VarType var, DirType dir, double *a, double *b, double *c, double *d, int n, TimeLayer3D *cur, TimeLayer3D *temp)
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
				a[p] = - temp->U->elem(i+p, j, k) / (2 * grid->dx) - vis_dx2; 
				b[p] = 3 / dt  +  2 * vis_dx2; 
				c[p] = temp->U->elem(i+p, j, k) / (2 * grid->dx) - vis_dx2; 
				
				switch (var)	
				{
				case type_U: d[p] = cur->U->elem(i+p, j, k) * 3 / dt - params.v_T * temp->T->d_x(i+p, j, k); break;
				case type_V: d[p] = cur->V->elem(i+p, j, k) * 3 / dt; break;
				case type_W: d[p] = cur->W->elem(i+p, j, k) * 3 / dt; break;
				case type_T: d[p] = cur->T->elem(i+p, j, k) * 3 / dt + params.t_phi * temp->DissFuncX(i+p, j, k); break;
				}	
				break;

			case Y:
				a[p] = - temp->V->elem(i, j+p, k) / (2 * grid->dy) - vis_dy2; 
				b[p] = 3 / dt  +  2 * vis_dy2; 
				c[p] = temp->V->elem(i, j+p, k) / (2 * grid->dy) - vis_dy2; 
				
				switch (var)	
				{
				case type_U: d[p] = cur->U->elem(i, j+p, k) * 3 / dt; break;
				case type_V: d[p] = cur->V->elem(i, j+p, k) * 3 / dt - params.v_T * temp->T->d_y(i, j+p, k); break;
				case type_W: d[p] = cur->W->elem(i, j+p, k) * 3 / dt; break;
				case type_T: d[p] = cur->T->elem(i, j+p, k) * 3 / dt + params.t_phi * temp->DissFuncY(i, j+p, k); break;
				}
				break;

			case Z:
				a[p] = - temp->W->elem(i, j, k+p) / (2 * grid->dz) - vis_dz2; 
				b[p] = 3 / dt  +  2 * vis_dz2; 
				c[p] = temp->W->elem(i, j, k+p) / (2 * grid->dz) - vis_dz2; 
				
				switch (var)	
				{
				case type_U: d[p] = cur->U->elem(i, j, k+p) * 3 / dt; break;
				case type_V: d[p] = cur->V->elem(i, j, k+p) * 3 / dt; break;
				case type_W: d[p] = cur->W->elem(i, j, k+p) * 3 / dt - params.v_T * temp->T->d_z(i, j, k+p); break;
				case type_T: d[p] = cur->T->elem(i, j, k+p) * 3 / dt + params.t_phi * temp->DissFuncZ(i, j, k+p); break;
				}
				break;
			}
		}
	}

	void AdiSolver3D::ApplyBC0(int i, int j, int k, VarType var, double &b0, double &c0, double &d0)
	{
		if ((var == type_T && grid->GetBC_temp(i, j, k) == BC_FREE) ||
			(var != type_T && grid->GetBC_vel(i, j, k) == BC_FREE))
		{
			// free
			b0 = 1.0; 
			c0 = -1.0; 
			d0 = 0.0; 
		}
		else
		{
			// no-slip
			b0 = 1.0; 
			c0 = 0.0; 
			switch (var)
			{
			case type_U: d0 = grid->GetVel(i, j, k).x; break;
			case type_V: d0 = grid->GetVel(i, j, k).y; break;
			case type_W: d0 = grid->GetVel(i, j, k).z; break;
			case type_T: d0 = grid->GetT(i, j, k); break;
			}
		}
	}

	void AdiSolver3D::ApplyBC1(int i, int j, int k, VarType var, double &a1, double &b1, double &d1)
	{
		if ((var == type_T && grid->GetBC_temp(i, j, k) == BC_FREE) ||
			(var != type_T && grid->GetBC_vel(i, j, k) == BC_FREE))
		{
			// free
			a1 = 1.0; 
			b1 = -1.0; 
			d1 = 0.0;
		}
		else
		{
			// no-slip
			a1 = 0.0; 
			b1 = 1.0; 
			switch (var)
			{
			case type_U: d1 = grid->GetVel(i, j, k).x; break;
			case type_V: d1 = grid->GetVel(i, j, k).y; break;
			case type_W: d1 = grid->GetVel(i, j, k).z; break;
			case type_T: d1 = grid->GetT(i, j, k); break;
			}
		}
	}	
}