#pragma once

#include <stdio.h>

#define INF				1e10

namespace FluidSolver
{	
	enum RetStatus { OK, ERR };

	struct Vec2D
	{
		double x, y; 

		Vec2D() : x(0.0), y(0.0) { }
		Vec2D(double _x, double _y) : x(_x), y(_y) { }
		Vec2D(Vec2D &vec) : x(vec.x), y(vec.y) { }
	};

	struct Point2D 
	{ 
		double x, y; 

		Point2D() : x(0.0), y(0.0) { }
		Point2D(double _x, double _y) : x(_x), y(_y) { }
	};

	struct BBox2D
	{
		Point2D pMin, pMax;

		BBox2D() { Clear(); }
		
		void AddPoint(Point2D p)
		{
			if (p.x < pMin.x) pMin.x = p.x;
			if (p.y < pMin.y) pMin.y = p.y;
			if (p.x > pMax.x) pMax.x = p.x;
			if (p.y > pMax.y) pMax.y = p.y;
		}

		void Clear()
		{
			pMin.x = pMin.y = INF; 
			pMax.x = pMax.y = -INF;
		}
	};

	struct FluidParams
	{
		double Re, Pr, lambda;

		FluidParams() { }
		FluidParams(double _Re, double _Pr, double _lambda) : Re(_Re), Pr(_Pr), lambda(_lambda) { }
	};

	static void OutputResultHeader(FILE *file, BBox2D *bbox, int outdimx, int outdimy)
	{
		fprintf(file, "%.2f %.2f %.2f %.2f\n", bbox->pMin.x * 1000, bbox->pMin.y * 1000, bbox->pMax.x * 1000, bbox->pMax.y * 1000);

		float ddx = (float)(bbox->pMax.x - bbox->pMin.x) / outdimx;
		float ddy = (float)(bbox->pMax.y - bbox->pMin.y) / outdimy;
		fprintf(file, "%.2f %.2f %i %i\n", ddx * 1000, ddy * 1000, outdimx, outdimy);
	}

	static void OutputResult(FILE* file, Vec2D *v, double *T, int dimx, int dimy, float timeValue)
	{
		fprintf(file, "%.3f\n", timeValue);
		for (int j = 0; j < dimy; j++)
		{
			for (int i = 0; i < dimx; i++)
				fprintf(file, "%.2f %.2f ", v[i * dimy + j].x * 10, v[i * dimy + j].y * 10);
			fprintf(file, "\n");
		}
	}

	static int LoadLastLayer(char *fileName, Vec2D *v, double *T, int dimx, int dimy, int frames)
	{
		FILE *file = NULL;
		if(!fopen_s(&file, fileName, "r"))
		{
			int indimx, indimy, frame;
			fscanf_s(file, "%i %i %i", &frame, &indimx, &indimy);
			if (indimx != dimx || indimy != dimy || frame <= 0 || frame > frames) 
			{
				fclose(file);
				return 0;
			}
			
			for (int j = 0; j < dimy; j++)
				for (int i = 0; i < dimx; i++)
				{
					float vx, vy, t;
					fscanf_s(file, "%f %f %f", &vx, &vy, &t);
					v[i * dimy + j].x = (double)vx;
					v[i * dimy + j].y = (double)vy;
					T[i * dimy + j] = (double)t;
				}
			
			fclose(file);
			return frame;
		}
		else 
			return 0;
	}

	static void SaveLastLayer(char *fileName, int frame, Vec2D *v, double *T, int dimx, int dimy)
	{
		FILE *file = NULL;
		fopen_s(&file, fileName, "w");
		fprintf(file, "%i\n", frame);	
		fprintf(file, "%i %i\n", dimx, dimy);
		for (int j = 0; j < dimy; j++)
		{
			for (int i = 0; i < dimx; i++)
				fprintf(file, "%f %f %f ", v[i * dimy + j].x, v[i * dimy + j].y, T[i * dimy + j]);
			fprintf(file, "\n");
		}
		fclose(file);
	}

	static void PrintTimeStepInfo(int i, int j, int frames, int subframes, float elapsed_time)
	{
		const int percent = frames * subframes;
		int step = (j + i * subframes) * 100;
		
		float time_left_sec = ((frames-i-1) * subframes + subframes-j-1) * elapsed_time;
		int time_h = ((int)time_left_sec) / 3600;
		int time_m = (((int)time_left_sec) / 60) % 60;
		int time_s = ((int)time_left_sec) % 60;
		
		printf(" frame %i\tsubstep %i\t%i%%\t(%i h %i m %i s left)\n", i, j, step / percent, time_h, time_m, time_s);
	}
}
