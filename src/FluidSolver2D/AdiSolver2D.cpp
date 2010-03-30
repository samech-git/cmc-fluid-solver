#include "AdiSolver2D.h"

namespace FluidSolver2D
{
	void AdiSolver2D::Init(Grid2D* _grid, FluidParams &_params)
	{
		grid = _grid;

		dimx = grid->dimx;
		dimy = grid->dimy;

		params = _params;

		cur = new TimeLayer2D(grid->dimx, grid->dimy, grid->dx, grid->dy);
		half = new TimeLayer2D(grid->dimx, grid->dimy, grid->dx, grid->dy);
		next = new TimeLayer2D(grid->dimx, grid->dimy, grid->dx, grid->dy);

		temp = new TimeLayer2D(grid->dimx, grid->dimy, grid->dx, grid->dy);
		next_local = new TimeLayer2D(grid->dimx, grid->dimy, grid->dx, grid->dy);

		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				switch(grid->GetType(i, j))
				{
				case IN:
				case BOUND:
				case VALVE:
				case OUT:
					cur->U(i, j) = grid->GetData(i, j).vel.x;
					cur->V(i, j) = grid->GetData(i, j).vel.y;
					cur->T(i, j) = grid->GetData(i, j).T;
					break;
				}
	}

	void AdiSolver2D::UpdateSegment(double *x, Segment2D seg, VarType var, TimeLayer2D *layer)
	{
		int i = seg.posx;
		int j = seg.posy;

		for (int t = 0; t < seg.size; t++)
		{
			switch (var)
			{
			case type_U: layer->U(i, j) = x[t]; break;
			case type_V: layer->V(i, j) = x[t]; break;
			case type_T: layer->T(i, j) = x[t]; break;
			}

			switch (seg.dir)
			{
			case X: i++; break;
			case Y: j++; break;
			}
		}
	}

	void AdiSolver2D::ApplyBC0(int i, int j, VarType var, double &b0, double &c0, double &d0)
	{
		switch (grid->GetData(i, j).type)
		{
		case NONE: 
		case NOSLIP: 
			b0 = 1.0; 
			c0 = 0.0; 
			switch (var)
			{
			case type_U: d0 = grid->GetData(i, j).vel.x; break;
			case type_V: d0 = grid->GetData(i, j).vel.y; break;
			case type_T: d0 = grid->GetData(i, j).T; break;
			}
			break;
		case FREE: 
			b0 = 1.0; 
			c0 = -1.0; 
			d0 = 0.0; 
			break;
		}
	}

	void AdiSolver2D::ApplyBC1(int i, int j, VarType var, double &a1, double &b1, double &d1)
	{
		switch (grid->GetData(i, j).type)
		{
		case NONE:
		case NOSLIP: 
			a1 = 0.0; 
			b1 = 1.0; 
			switch (var)
			{
			case type_U: d1 = grid->GetData(i, j).vel.x; break;
			case type_V: d1 = grid->GetData(i, j).vel.y; break;
			case type_T: d1 = grid->GetData(i, j).T; break;
			}
			break;
		case FREE: 
			a1 = 1.0; 
			b1 = -1.0; 
			d1 = 0.0;
			break;
		}
	}

	void AdiSolver2D::BuildMatrix(double dt, int i, int j, VarType var, DirType dir, double *a, double *b, double *c, double *d, int n, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *temp_local)
	{
		double v_vis_dx2 = params.v_vis / (grid->dx * grid->dx);
		double t_vis_dx2 = params.t_vis / (grid->dx * grid->dx);

		double v_vis_dy2 = params.v_vis / (grid->dy * grid->dy);
		double t_vis_dy2 = params.t_vis / (grid->dy * grid->dy);

		for (int p = 1; p < n-1; p++)
		{
			switch (dir)
			{
			case X:
				switch (var)	
				{
				case type_U:
					a[p] = - temp_local->U(i+p, j) / (2 * grid->dx) - v_vis_dx2; 
					b[p] = 1 / dt  +  2 * v_vis_dx2; 
					c[p] = temp_local->U(i+p, j) / (2 * grid->dx) - v_vis_dx2; 
					d[p] = cur->U(i+p, j) / dt - params.v_T * temp_local->Tx(i+p, j);
					break;
				case type_V:
					a[p] = - temp_local->U(i+p, j) / (2 * grid->dx) - v_vis_dx2; 
					b[p] = 1 / dt  +  2 * v_vis_dx2; 
					c[p] = temp_local->U(i+p, j) / (2 * grid->dx) - v_vis_dx2; 
					d[p] = cur->V(i+p, j) / dt;
					break;
				case type_T:
					a[p] = - temp_local->U(i+p, j) / (2 * grid->dx) - t_vis_dx2; 
					b[p] = 1 / dt  +  2 * t_vis_dx2; 
					c[p] = temp_local->U(i+p, j) / (2 * grid->dx) - t_vis_dx2; 
					d[p] = cur->T(i+p, j) / dt + params.t_phi * temp_local->DissFuncX(i+p, j);
					break;
				}
				break;
			case Y:
				switch (var)	
				{
				case type_U:
					a[p] = - temp_local->V(i, j+p) / (2 * grid->dy) - v_vis_dy2; 
					b[p] = 1 / dt  +  2 * v_vis_dy2; 
					c[p] = temp_local->V(i, j+p) / (2 * grid->dy) - v_vis_dy2; 
					d[p] = cur->U(i, j+p) / dt;
					break;
				case type_V:
					a[p] = - temp_local->V(i, j+p) / (2 * grid->dy) - v_vis_dy2; 
					b[p] = 1 / dt  +  2 * v_vis_dy2; 
					c[p] = temp_local->V(i, j+p) / (2 * grid->dy) - v_vis_dy2; 
					d[p] = cur->V(i, j+p) / dt - params.v_T * temp_local->Ty(i, j+p);
					break;
				case type_T:
					a[p] = - temp_local->V(i, j+p) / (2 * grid->dy) - t_vis_dy2; 
					b[p] = 1 / dt  +  2 * t_vis_dy2; 
					c[p] = temp_local->V(i, j+p) / (2 * grid->dy) - t_vis_dy2; 
					d[p] = cur->T(i, j+p) / dt + params.t_phi * temp_local->DissFuncY(i, j+p);
					break;
				}
				break;
			}
		}
	}

	void AdiSolver2D::SolveSegment(double dt, Segment2D seg, VarType var, DirType dir, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *temp_local, TimeLayer2D *next_local)
	{		
		int n = seg.size;

		a = new double[n];
		b = new double[n];
		c = new double[n];
		d = new double[n];
		x = new double[n];

		ApplyBC0(seg.posx, seg.posy, var, b[0], c[0], d[0]);
		BuildMatrix(dt, seg.posx, seg.posy, var, dir, a, b, c, d, n, cur, temp, temp_local);
		ApplyBC1(seg.endx, seg.endy, var, a[n-1], b[n-1], d[n-1]);
		SolveTridiagonal(a, b, c, d, x, n);
		UpdateSegment(x, seg, var, next_local);

		delete [] a;
		delete [] b;
		delete [] c;
		delete [] d;
		delete [] x;
	}

	void AdiSolver2D::SolveDirection(double dt, int num_local, vector<Segment2D> &list, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next)
	{
		TimeLayer2D *temp_local = new TimeLayer2D(grid->dimx, grid->dimy, grid->dx, grid->dy);
		temp->CopyAllto(grid, temp_local);

		DirType dir = list[0].dir;
		for (int it = 0; it < num_local; it++)
		{
			for (size_t s = 0; s < list.size(); s++)
			{		
				SolveSegment(dt, list[s], type_U, dir, cur, temp, temp_local, next_local);
				SolveSegment(dt, list[s], type_V, dir, cur, temp, temp_local, next_local);
				SolveSegment(dt, list[s], type_T, dir, cur, temp, temp_local, next_local);
			}

			// update non-linear
			if (it == 0) next_local->CopyAllto(grid, temp_local, IN);
				else next_local->MergeAllto(grid, temp_local, IN);
		}

		temp_local->CopyAllto(grid, temp, IN);
		next_local->CopyAllto(grid, next, IN);
		delete temp_local;
	}

	void AdiSolver2D::CreateSegments()
	{
		listX.clear();
		for (int i = 0; i < dimx; i++)
		{
			Segment2D seg;
			seg.posx = i;
			seg.dir = Y;

			int j = 0;
			while (j < dimy && grid->GetType(i, j) == OUT) j++;
			while (j+1 < dimy && grid->GetType(i, j+1) != IN) j++;
			if (j+1 >= dimy) continue;
			seg.posy = j;

			j = dimy-1;
			while (j >= 0 && grid->GetType(i, j) == OUT) j--;
			while (j-1 >= 0 && grid->GetType(i, j-1) != IN) j--;
			seg.size = j - seg.posy + 1;
			
			seg.endx = i;
			seg.endy = j;

			listX.push_back(seg);
		}

		listY.clear();
		for (int j = 0; j < dimy; j++)
		{
			Segment2D seg;
			seg.posy = j;
			seg.dir = X;

			int i = 0;
			while (i < dimx && grid->GetType(i, j) == OUT) i++;
			while (i+1 < dimx && grid->GetType(i+1, j) != IN) i++;
			if (i+1 >= dimx) continue;
			seg.posx = i;

			i = dimx-1;
			while (i >= 0 && grid->GetType(i, j) == OUT) i--;
			while (i-1 >= 0 && grid->GetType(i-1, j) != IN) i--;
			seg.size = i - seg.posx + 1;

			seg.endx = i;
			seg.endy = j;

			listY.push_back(seg);
		}
	}

	void AdiSolver2D::TimeStep(double dt, int num_global, int num_local)
	{
		CreateSegments();

		cur->CopyAllto(grid, next);
		cur->CopyAllto(grid, half);
		cur->CopyAllto(grid, temp);

		// do global iterations
		int it;
		double err = next->EvalDivError(grid);

		for (it = 0; (it < num_global) || (err > ERR_THRESHOLD); it++)
		{
			// alternating directions
			SolveDirection(dt, num_local, listY, cur, temp, half);
			SolveDirection(dt, num_local, listX, half, temp, next);

			err = next->EvalDivError(grid);
			
			// update non-linear parameters
			if (it == 0) next->CopyAllto(grid, temp, IN);
				else next->MergeAllto(grid, temp, IN);
			
			if (it > MAX_GLOBAL_ITERS) 
			{
				printf("\nExceeded max number of iterations (%i)\n", MAX_GLOBAL_ITERS); 
				exit(1); 
			}

			if (err > ERR_THRESHOLD * 10)
			{
				printf("\nError is too big!\n", err);
				exit(1);
			}
		}

		ClearOutterCells();

		// output error
		printf("\rerr = %.4f,", err);

		// copy result to current layer
		next->CopyAllto(grid, cur);
	}

	void AdiSolver2D::FreeMemory()
	{
		if (cur != NULL) delete cur;
		if (temp != NULL) delete temp;
		if (half != NULL) delete half;
		if (next_local != NULL) delete next_local;
		if (next != NULL) delete next;
	}

	AdiSolver2D::AdiSolver2D()
	{
		grid = NULL;

		cur = NULL;
		temp = NULL;
		half = NULL;
		next_local = NULL;
		next = NULL;
	}

	AdiSolver2D::~AdiSolver2D()
	{
		FreeMemory();
	}
}