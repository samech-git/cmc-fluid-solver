#pragma once

#include "Grid2D.h"

namespace FluidSolver2D
{
	struct TimeLayer2D
	{
		int dimx, dimy;
		double dx, dy;

		double& U(int i, int j)
		{
			return u[i * dimy + j];
		}

		double& V(int i, int j)
		{
			return v[i * dimy + j];
		}

		double& T(int i, int j)
		{
			return t[i * dimy + j];
		}

		// U derivatives
		inline double Ux(int i, int j)	{ return (U(i+1, j) - U(i-1, j)) / (2 * dx); }
		inline double Uy(int i, int j)	{ return (U(i, j+1) - U(i, j-1)) / (2 * dy); }
		inline double Uxx(int i, int j) { return (U(i+1, j) - 2 * U(i, j) + U(i-1, j)) / (dx * dx); }
		inline double Uyy(int i, int j) { return (U(i, j+1) - 2 * U(i, j) + U(i, j-1)) / (dy * dy); }
		
		// V derivatives
		inline double Vx(int i, int j)	{ return (V(i+1, j) - V(i-1, j)) / (2 * dx); }
		inline double Vy(int i, int j)	{ return (V(i, j+1) - V(i, j-1)) / (2 * dy); }
		inline double Vxx(int i, int j) { return (V(i+1, j) - 2 * V(i, j) + V(i-1, j)) / (dx * dx); }
		inline double Vyy(int i, int j) { return (V(i, j+1) - 2 * V(i, j) + V(i, j-1)) / (dy * dy); }

		// T derivatives
		inline double Tx(int i, int j)	{ return (T(i+1, j) - T(i-1, j)) / (2 * dx); }
		inline double Ty(int i, int j)	{ return (T(i, j+1) - T(i, j-1)) / (2 * dy); }
		inline double Txx(int i, int j) { return (T(i+1, j) - 2 * T(i, j) + T(i-1, j)) / (dx * dx); }
		inline double Tyy(int i, int j) { return (T(i, j+1) - 2 * T(i, j) + T(i, j-1)) / (dy * dy); }
		
		// other functions

		double DissFuncX(int i, int j)
		{
			double ux = Ux(i, j);
			double vx = Vx(i, j);
			double uy = Uy(i, j);
			double vy = Vy(i, j);

			return 2 * ux * ux + vx * vx + uy * vx;
		}

		double DissFuncY(int i, int j)
		{
			double ux = Ux(i, j);
			double vx = Vx(i, j);
			double uy = Uy(i, j);
			double vy = Vy(i, j);

			return uy * uy + 2 * vy * vy + vx * uy;
		}

		double DissFunc(int i, int j)
		{
			return DissFuncX(i, j) + DissFuncY(i, j);
		}

		double EvalDivError(Grid2D *grid)
		{
			double err = 0.0;
			int count = 0;
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					if (grid->GetType(i, j) == IN && grid->GetType(i+1, j) == IN && grid->GetType(i, j+1) == IN && grid->GetType(i+1, j+1) == IN)
					{
						double tx = dy * (U(i+1, j) - U(i, j)) + (U(i+1, j+1) - U(i, j+1)) / 2;
						double ty = dx * (V(i, j+1) - V(i, j)) + (V(i+1, j+1) - V(i+1, j)) / 2;
						err += abs(tx + ty);
						count++;
					}
			return err / count;
		}

		void CopyUto(Grid2D *grid, TimeLayer2D *dest, CellType type)
		{
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					if (grid->GetType(i, j) == type)
						dest->U(i, j) = U(i, j);
		}

		void MergeUto(Grid2D *grid, TimeLayer2D *dest, CellType type)
		{
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					if (grid->GetType(i, j) == type)
						dest->U(i, j) = (dest->U(i, j) + U(i, j)) / 2;
		}

		void CopyVto(Grid2D *grid, TimeLayer2D *dest, CellType type)
		{
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					if (grid->GetType(i, j) == type)
						dest->V(i, j) = V(i, j);
		}

		void MergeVto(Grid2D *grid, TimeLayer2D *dest, CellType type)
		{
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					if (grid->GetType(i, j) == type)
						dest->V(i, j) = (dest->V(i, j) + V(i, j)) / 2;
		}

		void CopyTto(Grid2D *grid, TimeLayer2D *dest, CellType type)
		{
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					if (grid->GetType(i, j) == type)
						dest->T(i, j) = T(i, j);
		}

		void MergeTto(Grid2D *grid, TimeLayer2D *dest, CellType type)
		{
			for (int i = 0; i < dimx-1; i++)
				for (int j = 0; j < dimy-1; j++)
					if (grid->GetType(i, j) == type)
						dest->T(i, j) = (dest->T(i, j) + T(i, j)) / 2;
		}

		void MergeAllto(Grid2D *grid, TimeLayer2D *dest, CellType type)
		{
			MergeUto(grid, dest, type);
			MergeVto(grid, dest, type);
			MergeTto(grid, dest, type);
		}

		void CopyAllto(Grid2D *grid, TimeLayer2D *dest)
		{
			CopyAllto(grid, dest, IN);
			CopyAllto(grid, dest, OUT);
			CopyAllto(grid, dest, BOUND);
			CopyAllto(grid, dest, VALVE);
		}

		void CopyAllto(Grid2D *grid, TimeLayer2D *dest, CellType type)
		{
			CopyUto(grid, dest, type);
			CopyVto(grid, dest, type);
			CopyTto(grid, dest, type);
		}

		TimeLayer2D(int _dimx, int _dimy, double _dx, double _dy) : dimx(_dimx), dimy(_dimy), dx(_dx), dy(_dy)
		{
			u = new double[dimx * dimy];
			v = new double[dimx * dimy];
			t = new double[dimx * dimy];
		}
		
		~TimeLayer2D()
		{
			delete [] u;
			delete [] v;
			delete [] t;
		}

	private:	
		double *u, *v, *t;
	};
}