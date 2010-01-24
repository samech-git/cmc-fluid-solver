#pragma once

#include <stdio.h>

namespace FluidSolver
{
	enum RetStatus { OK, ERR };

	struct Vec2D
	{
		double x, y; 

		Vec2D() : x(0.0), y(0.0) { }
		Vec2D(double _x, double _y) : x(_x), y(_y) { }
	};

	struct FluidParams
	{
		double Re, Pr, lambda;

		FluidParams() { }
		FluidParams(double _Re, double _Pr, double _lambda) : Re(_Re), Pr(_Pr), lambda(_lambda) { }
	};

	static void TestPrintResult(int dimx, int dimy, Vec2D *v, double *T)
	{
		printf("writing results..\n");
		FILE *file = NULL;
		fopen_s(&file, "results.txt", "w");
		fprintf(file, "u:\n");
		for (int i = 0; i < dimx; i++)
		{
			for (int j = 0; j < dimy; j++)
				fprintf(file, "%.2f ", v[i * dimy + j].x);
			fprintf(file, "\n");
		}
		fprintf(file, "v:\n");
		for (int i = 0; i < dimx; i++)
		{
			for (int j = 0; j < dimy; j++)
				fprintf(file, "%.2f ", v[i * dimy + j].y);
			fprintf(file, "\n");
		}
		fprintf(file, "T:\n");
		for (int i = 0; i < dimx; i++)
		{
			for (int j = 0; j < dimy; j++)
				fprintf(file, "%.2f ", T[i * dimy + j]);
			fprintf(file, "\n");
		}
	}


	static void ShiferTestPrintResult(FILE* file, Vec2D *v, double *T, int dimx, int dimy)
	{
		for (int j = 0; j < dimy; j++)
		{
			for (int i = 0; i < dimx; i++)
				fprintf(file, "%.2f %.2f ", v[i * dimy + j].x, v[i * dimy + j].y);
			fprintf(file, "\n");
		}
	}

	static int LoadLastLayer(char *fileName, Vec2D **v, double **T, int dimx, int dimy)
	{
		*v = new Vec2D[dimx * dimy];
		*T = new double[dimx * dimy];

		return 0;
	}

	static void PrintLast(FILE* file, Vec2D *v, double *T, int dimx, int dimy)
	{
		for (int j = 0; j < dimy; j++)
		{
			for (int i = 0; i < dimx; i++)
				fprintf(file, "%f %f %f", v[i * dimy + j].x, v[i * dimy + j].y, T[i * dimy + j]);
			fprintf(file, "\n");
		}
	}
}
