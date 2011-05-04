/*
 *  Copyright 2010-2011 Nikolai Sakharnykh
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#pragma once

#include "Grid2D.h"

namespace FluidSolver2D
{
	struct ScalarField2D
	{
		int dimx, dimy;
		FTYPE dx, dy;
		
		FTYPE& U(int i, int j)
		{
			return u[i * dimy + j];
		}

		// U derivatives
		inline FTYPE Ux(int i, int j)	{ return (U(i+1, j) - U(i-1, j)) / (2 * dx); }
		inline FTYPE Uy(int i, int j)	{ return (U(i, j+1) - U(i, j-1)) / (2 * dy); }
		inline FTYPE Uxx(int i, int j) { return (U(i+1, j) - 2 * U(i, j) + U(i-1, j)) / (dx * dx); }
		inline FTYPE Uyy(int i, int j) { return (U(i, j+1) - 2 * U(i, j) + U(i, j-1)) / (dy * dy); }
		
		inline Vec2D grad(int i, int j) { return Vec2D(Ux(i, j), Uy(i, j)); }

		ScalarField2D(int _dimx, int _dimy, FTYPE _dx, FTYPE _dy) : dimx(_dimx), dimy(_dimy), dx(_dx), dy(_dy)
		{
			u = new FTYPE[dimx * dimy];
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
		FTYPE *u;
	};
}