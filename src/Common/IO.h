#pragma once

#include <stdio.h>

#define MAX_STR_SIZE	255

namespace Common
{	
	static void OutputResultHeader(FILE *file, BBox2D *bbox, int outdimx, int outdimy)
	{
		fprintf(file, "%.2f %.2f %.2f %.2f\n", bbox->pMin.x * 1000, bbox->pMin.y * 1000, bbox->pMax.x * 1000, bbox->pMax.y * 1000);

		float ddx = (float)(bbox->pMax.x - bbox->pMin.x) / outdimx;
		float ddy = (float)(bbox->pMax.y - bbox->pMin.y) / outdimy;
		fprintf(file, "%.2f %.2f %i %i\n", ddx * 1000, ddy * 1000, outdimx, outdimy);
	}

	static void OutputResult(FILE* file, Vec2D *v, double *T, int dimx, int dimy, float timeValue)
	{
		fprintf(file, "%.5f\n", timeValue);
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

	static void PrintTimeStepInfo(int frame, int subframe, double cur_time, double max_time, float elapsed_time)
	{
		//printf("elapsed sec: %.4f\n", elapsed_time);
		//const int percent = frames * subframes * cycles;
		//int step = (j + i * subframes) * 100;
		//int perres = step / percent;
		float perresf = (float)cur_time * 100 / (float)max_time;

		if (perresf < 2)
		{
			printf(" frame %i\tsubstep %i\t%i%%\t(----- left)", frame, subframe, (int)perresf);
		}
		else
		{
			//float time_left_sec = ((frames-i-1) * subframes + subframes-j-1) * elapsed_time;
			float time_left_sec = elapsed_time * (100 - perresf) / perresf;
			int time_h = ((int)time_left_sec) / 3600;
			int time_m = (((int)time_left_sec) / 60) % 60;
			int time_s = ((int)time_left_sec) % 60;

			printf(" frame %i\tsubstep %i\t%i%%\t(%i h %i m %i s left)", frame, subframe, (int)perresf, time_h, time_m, time_s);
		}
	}

	static void FindFile(char *path, char *filename, bool checkExist = true)
	{
		FILE *file = NULL;

		sprintf_s(path, MAX_STR_SIZE, "%s", filename);
		fopen_s(&file, path, "r");
		if (!file && checkExist) 
		{
			sprintf_s(path, MAX_STR_SIZE, "..\\..\\data\\%s", filename);
			fopen_s(&file, path, "r");
			if (!file) { printf("cannot find the file: \"%s\"\n", filename); }
				else fclose(file);	
		}
		else
			if (checkExist) fclose(file);
	}
}
