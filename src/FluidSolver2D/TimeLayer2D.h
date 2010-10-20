#pragma once

#include "Grid2D.h"

namespace FluidSolver2D
{
	struct TimeLayer2D
	{
		int dimx, dimy;
		FTYPE dx, dy;

		FTYPE& U(int i, int j)
		{
			return u[i * dimy + j];
		}

		FTYPE& V(int i, int j)
		{
			return v[i * dimy + j];
		}

		FTYPE& T(int i, int j)
		{
			return t[i * dimy + j];
		}

		// U derivatives
		inline FTYPE Ux(int i, int j)	{ return (U(i+1, j) - U(i-1, j)) / (2 * dx); }
		inline FTYPE Uy(int i, int j)	{ return (U(i, j+1) - U(i, j-1)) / (2 * dy); }
		inline FTYPE Uxx(int i, int j) { return (U(i+1, j) - 2 * U(i, j) + U(i-1, j)) / (dx * dx); }
		inline FTYPE Uyy(int i, int j) { return (U(i, j+1) - 2 * U(i, j) + U(i, j-1)) / (dy * dy); }
		
		// V derivatives
		inline FTYPE Vx(int i, int j)	{ return (V(i+1, j) - V(i-1, j)) / (2 * dx); }
		inline FTYPE Vy(int i, int j)	{ return (V(i, j+1) - V(i, j-1)) / (2 * dy); }
		inline FTYPE Vxx(int i, int j) { return (V(i+1, j) - 2 * V(i, j) + V(i-1, j)) / (dx * dx); }
		inline FTYPE Vyy(int i, int j) { return (V(i, j+1) - 2 * V(i, j) + V(i, j-1)) / (dy * dy); }

		// T derivatives
		inline FTYPE Tx(int i, int j)	{ return (T(i+1, j) - T(i-1, j)) / (2 * dx); }
		inline FTYPE Ty(int i, int j)	{ return (T(i, j+1) - T(i, j-1)) / (2 * dy); }
		inline FTYPE Txx(int i, int j) { return (T(i+1, j) - 2 * T(i, j) + T(i-1, j)) / (dx * dx); }
		inline FTYPE Tyy(int i, int j) { return (T(i, j+1) - 2 * T(i, j) + T(i, j-1)) / (dy * dy); }
		
		// other functions

		FTYPE DissFuncX(int i, int j)
		{
			FTYPE ux = Ux(i, j);
			FTYPE vx = Vx(i, j);
			FTYPE uy = Uy(i, j);
			FTYPE vy = Vy(i, j);

			return 2 * ux * ux + vx * vx + uy * vx;
		}

		FTYPE DissFuncY(int i, int j)
		{
			FTYPE ux = Ux(i, j);
			FTYPE vx = Vx(i, j);
			FTYPE uy = Uy(i, j);
			FTYPE vy = Vy(i, j);

			return uy * uy + 2 * vy * vy + vx * uy;
		}

		FTYPE DissFunc(int i, int j)
		{
			return DissFuncX(i, j) + DissFuncY(i, j);
		}

		double EvalDivError(Grid2D *grid)
		{
			FTYPE err = 0.0;
			int count = 0;
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					if (grid->GetType(i, j) == NODE_IN && grid->GetType(i+1, j) == NODE_IN && grid->GetType(i, j+1) == NODE_IN && grid->GetType(i+1, j+1) == NODE_IN)
					{
						FTYPE tx = dy * (U(i+1, j) - U(i, j)) + (U(i+1, j+1) - U(i, j+1)) / 2;
						FTYPE ty = dx * (V(i, j+1) - V(i, j)) + (V(i+1, j+1) - V(i+1, j)) / 2;
						err += abs(tx + ty);
						count++;
					}
			return err / count;
		}

		void CopyUto(Grid2D *grid, TimeLayer2D *dest, NodeType type)
		{
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					if (grid->GetType(i, j) == type)
						dest->U(i, j) = U(i, j);
		}

		void MergeUto(Grid2D *grid, TimeLayer2D *dest, NodeType type)
		{
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					if (grid->GetType(i, j) == type)
						dest->U(i, j) = (dest->U(i, j) + U(i, j)) / 2;
		}

		void CopyVto(Grid2D *grid, TimeLayer2D *dest, NodeType type)
		{
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					if (grid->GetType(i, j) == type)
						dest->V(i, j) = V(i, j);
		}

		void MergeVto(Grid2D *grid, TimeLayer2D *dest, NodeType type)
		{
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					if (grid->GetType(i, j) == type)
						dest->V(i, j) = (dest->V(i, j) + V(i, j)) / 2;
		}

		void CopyTto(Grid2D *grid, TimeLayer2D *dest, NodeType type)
		{
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					if (grid->GetType(i, j) == type)
						dest->T(i, j) = T(i, j);
		}

		void MergeTto(Grid2D *grid, TimeLayer2D *dest, NodeType type)
		{
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					if (grid->GetType(i, j) == type)
						dest->T(i, j) = (dest->T(i, j) + T(i, j)) / 2;
		}

		void MergeAllto(Grid2D *grid, TimeLayer2D *dest, NodeType type)
		{
			MergeUto(grid, dest, type);
			MergeVto(grid, dest, type);
			MergeTto(grid, dest, type);
		}

		void CopyAllto(Grid2D *grid, TimeLayer2D *dest)
		{
			CopyAllto(grid, dest, NODE_IN);
			CopyAllto(grid, dest, NODE_OUT);
			CopyAllto(grid, dest, NODE_BOUND);
			CopyAllto(grid, dest, NODE_VALVE);
		}

		void CopyAllto(Grid2D *grid, TimeLayer2D *dest, NodeType type)
		{
			CopyUto(grid, dest, type);
			CopyVto(grid, dest, type);
			CopyTto(grid, dest, type);
		}

		TimeLayer2D(int _dimx, int _dimy, FTYPE _dx, FTYPE _dy) : dimx(_dimx), dimy(_dimy), dx(_dx), dy(_dy)
		{
			u = new FTYPE[dimx * dimy];
			v = new FTYPE[dimx * dimy];
			t = new FTYPE[dimx * dimy];
		}
		
		~TimeLayer2D()
		{
			delete [] u;
			delete [] v;
			delete [] t;
		}

	private:	
		FTYPE *u, *v, *t;
	};
}