#pragma once

#include "Grid2D.h"

namespace FluidSolver2D
{
	struct ScalarField2D
	{
		int dimx, dimy;
		double dx, dy;
		
		double& U(int i, int j)
		{
			return u[i * dimy + j];
		}

		// U derivatives
		inline double Ux(int i, int j)	{ return (U(i+1, j) - U(i-1, j)) / (2 * dx); }
		inline double Uy(int i, int j)	{ return (U(i, j+1) - U(i, j-1)) / (2 * dy); }
		inline double Uxx(int i, int j) { return (U(i+1, j) - 2 * U(i, j) + U(i-1, j)) / (dx * dx); }
		inline double Uyy(int i, int j) { return (U(i, j+1) - 2 * U(i, j) + U(i, j-1)) / (dy * dy); }
		
		inline Vec2D grad(int i, int j) { return Vec2D(Ux(i, j), Uy(i, j)); }

		ScalarField2D(int _dimx, int _dimy, double _dx, double _dy) : dimx(_dimx), dimy(_dimy), dx(_dx), dy(_dy)
		{
			u = new double[dimx * dimy];
		}

		void ClearZero()
		{
			for (int i = 0; i < dimx * dimy; i++) 
				u[i] = 0.0;
		}
		
		~ScalarField2D()
		{
			delete [] u;
		}
	private:	
		double *u;
	};
}