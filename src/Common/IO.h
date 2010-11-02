#pragma once

#include <stdio.h>
#include <string.h>
#include <direct.h>

#include <string>
#include <iostream>
#include <sstream>

#define MAX_STR_SIZE	255

namespace Common
{	
	#pragma pack(push,2)
	typedef struct tagBitmapFileHeader {
        unsigned short		bfType;
        unsigned long		bfSize;
        unsigned short		bfReserved1;
        unsigned short		bfReserved2;
        unsigned long		bfOffBits;
	} BitmapFileHeader;
	#pragma pack(pop)

	typedef struct tagBitmapInfoHeader {
        unsigned long		biSize;
        long				biWidth;
        long				biHeight;
        unsigned short      biPlanes;
        unsigned short      biBitCount;
        unsigned long		biCompression;
        unsigned long		biSizeImage;
		long			    biXPelsPerMeter;
        long				biYPelsPerMeter;
        unsigned long		biClrUsed;
        unsigned long		biClrImportant;
	} BitmapInfoHeader;

	inline std::string stringify(int x)
	{
		std::ostringstream o;
		o << x;
		return o.str();
	}

	static void OutputResultHeader(const char *outputPath, BBox2D *bbox, int outdimx, int outdimy)
	{
		FILE *file = NULL;
		fopen_s(&file, outputPath, "w");

		fprintf(file, "%.2f %.2f %.2f %.2f\n", bbox->pMin.x * 1000, bbox->pMin.y * 1000, bbox->pMax.x * 1000, bbox->pMax.y * 1000);

		float ddx = (float)(bbox->pMax.x - bbox->pMin.x) / outdimx;
		float ddy = (float)(bbox->pMax.y - bbox->pMin.y) / outdimy;
		fprintf(file, "%.2f %.2f %i %i\n", ddx * 1000, ddy * 1000, outdimx, outdimy);
		
		fclose(file);
	}

	static void OutputResult(const char *outputPath, Vec2D *v, double *T, int dimx, int dimy, float timeValue)
	{
		FILE* file = NULL;
		fopen_s(&file, outputPath, "a");
		
		fprintf(file, "%.5f\n", timeValue);
		for (int j = 0; j < dimy; j++)
		{
			for (int i = 0; i < dimx; i++)
				fprintf(file, "%.2f %.2f ", v[i * dimy + j].x * 10, v[i * dimy + j].y * 10);
			fprintf(file, "\n");
		}

		fclose(file);
	}

	// output Z-slice, velocity projected onto XY
	static void OutputSliceResult(const char *outputPath, int z, Vec3D *v, double *T, int dimx, int dimy, int dimz, float timeValue)
	{
		FILE *file = NULL;
		fopen_s(&file, outputPath, "a");
			
		fprintf(file, "%.5f\n", timeValue);
		for (int j = 0; j < dimy; j++)
		{
			for (int i = 0; i < dimx; i++)
				fprintf(file, "%.2f %.2f ", v[i * dimy * dimz + j * dimz + z].x * 10, v[i * dimy * dimz + j * dimz + z].y * 10);
			fprintf(file, "\n");
		}

		fclose(file);
	}

	static void OutputNetCDFHeader3D(const char *outputPath, BBox2D *bbox, double depth, double timestep, double time, int outdimx, int outdimy, int outdimz)
	{
		FILE *file = NULL;
		fopen_s(&file, outputPath, "w");

		fprintf(file, "netcdf 3d_scalar_time_array {\n");
	
		// --- dimensions ---
		fprintf(file, "dimensions:\n");
		fprintf(file, "\tx = %i ;\n", outdimx);
		fprintf(file, "\ty = %i ;\n", outdimy);
		fprintf(file, "\tz = %i ;\n", outdimz);
		fprintf(file, "\ttime = UNLIMITED ;\n");
	
		// --- variables ---
		fprintf(file, "variables:\n");
		
		fprintf(file, "\tfloat x(x) ;\n");
		fprintf(file, "\t\tx:units = \"metres\" ;\n");
		fprintf(file, "\t\tx:actual_range = %.2ff, %.2ff ;\n", bbox->pMin.x, bbox->pMax.x);
		fprintf(file, "\t\tx:long_name = \"X coordinate\" ;\n");

		fprintf(file, "\tfloat y(y) ;\n");
		fprintf(file, "\t\ty:units = \"metres\" ;\n");
		fprintf(file, "\t\ty:actual_range = %.2ff, %.2ff ;\n", bbox->pMin.y, bbox->pMax.y);
		fprintf(file, "\t\ty:long_name = \"Y coordinate\" ;\n");

		fprintf(file, "\tfloat z(z) ;\n");
		fprintf(file, "\t\tz:units = \"metres\" ;\n");
		fprintf(file, "\t\tz:actual_range = %.2ff, 0.f ;\n", -depth);
		fprintf(file, "\t\tz:long_name = \"Z coordinate\" ;\n");

		fprintf(file, "\tdouble time(time) ;\n");
		fprintf(file, "\t\ttime:units = \"s\" ;\n");
		fprintf(file, "\t\ttime:actual_range = 0.f, %.2ff ;\n", time);
		fprintf(file, "\t\ttime:long_name = \"Time\" ;\n");

		fprintf(file, "\tdouble u(time, x, y, z) ;\n");
		fprintf(file, "\t\tu:units = \"m/s\" ;\n");
		fprintf(file, "\t\tu:actual_range = 0.f, 1.f ;\n");
		fprintf(file, "\t\tu:valid_range = 0.f, 1.f ;\n");
		fprintf(file, "\t\tu:long_name = \"U velocity\" ;\n");
		fprintf(file, "\t\tu:scale_factor =  1.f ;\n");
		fprintf(file, "\t\tu:var_desc = \"U velocity\",\n\t\t\t\"U\" ; \n");

		fprintf(file, "\t// global attributes\n");
		fprintf(file, "\t:Conventions = \"COARDS\" ;\n");
		fprintf(file, "\t:title = \"3D Time U velocity data from FluidSolver3D (http://code.google.com/p/cmc-fluid-solver/)\" ;\n");
		fprintf(file, "\t:history = \"created by using FluidSolver3D library\" ;\n");
		fprintf(file, "\t:description = \"Test data\" ;\n");
		fprintf(file, "\t:platform = \"Model\" ;\n");

		// --- data ---
		fprintf(file, "data:\n");

		float ddx = (float)(bbox->pMax.x - bbox->pMin.x) / outdimx;
		float ddy = (float)(bbox->pMax.y - bbox->pMin.y) / outdimy;
		float ddz = (float)depth / outdimz;

		fprintf(file, "x = ");
		for (int i = 0; i < outdimx-1; i++)
			fprintf(file, "%.2f, ", bbox->pMin.x + ddx * i);
		fprintf(file, "%.2f ;\n", bbox->pMin.x + ddx * outdimx);

		fprintf(file, "y = ");
		for (int i = 0; i < outdimy-1; i++)
			fprintf(file, "%.2f, ", bbox->pMin.y + ddy * i);
		fprintf(file, "%.2f ;\n", bbox->pMin.y + ddy * outdimy);

		fprintf(file, "z = ");
		for (int i = 0; i < outdimz-1; i++)
			fprintf(file, "%.2f, ", ddz * i);
		fprintf(file, "%.2f ;\n", ddz * outdimz);
		
		fprintf(file, "time = ");
		for (float cur = 0; cur < time; cur += (float)timestep)
			fprintf(file, "%.2f, ", cur);
		fprintf(file, "%.2f ;\n", time);

		fprintf(file, "u = \n");

		fclose(file);
	}

	static void OutputNetCDFHeader2D(const char *outputPath, BBox2D *bbox, double timestep, double time, int outdimx, int outdimy)
	{
		FILE *file = NULL;
		fopen_s(&file, outputPath, "w");

		fprintf(file, "netcdf 2d_scalar_time_array {\n");
	
		// --- dimensions ---
		fprintf(file, "dimensions:\n");
		fprintf(file, "\tx = %i ;\n", outdimx);
		fprintf(file, "\ty = %i ;\n", outdimy);
		fprintf(file, "\ttime = UNLIMITED ;\n");
	
		// --- variables ---
		fprintf(file, "variables:\n");
		
		fprintf(file, "\tfloat x(x) ;\n");
		fprintf(file, "\t\tx:units = \"metres\" ;\n");
		fprintf(file, "\t\tx:actual_range = %.2ff, %.2ff ;\n", bbox->pMin.x, bbox->pMax.x);
		fprintf(file, "\t\tx:long_name = \"X coordinate\" ;\n");

		fprintf(file, "\tfloat y(y) ;\n");
		fprintf(file, "\t\ty:units = \"metres\" ;\n");
		fprintf(file, "\t\ty:actual_range = %.2ff, %.2ff ;\n", bbox->pMin.y, bbox->pMax.y);
		fprintf(file, "\t\ty:long_name = \"Y coordinate\" ;\n");

		fprintf(file, "\tdouble time(time) ;\n");
		fprintf(file, "\t\ttime:units = \"s\" ;\n");
		fprintf(file, "\t\ttime:actual_range = 0.f, %.2ff ;\n", time);
		fprintf(file, "\t\ttime:long_name = \"Time\" ;\n");

		fprintf(file, "\tdouble u(time, x, y) ;\n");
		fprintf(file, "\t\tu:units = \"m/s\" ;\n");
		fprintf(file, "\t\tu:actual_range = 0.f, 1.f ;\n");
		fprintf(file, "\t\tu:valid_range = 0.f, 1.f ;\n");
		fprintf(file, "\t\tu:long_name = \"U velocity\" ;\n");
		fprintf(file, "\t\tu:scale_factor =  1.f ;\n");
		fprintf(file, "\t\tu:var_desc = \"U velocity\",\n\t\t\t\"U\" ; \n");

		fprintf(file, "\t// global attributes\n");
		fprintf(file, "\t:Conventions = \"COARDS\" ;\n");
		fprintf(file, "\t:title = \"2D Time U velocity data from FluidSolver2D (http://code.google.com/p/cmc-fluid-solver/)\" ;\n");
		fprintf(file, "\t:history = \"created by using FluidSolver2D library\" ;\n");
		fprintf(file, "\t:description = \"Test data\" ;\n");
		fprintf(file, "\t:platform = \"Model\" ;\n");

		// --- data ---
		fprintf(file, "data:\n");

		float ddx = (float)(bbox->pMax.x - bbox->pMin.x) / outdimx;
		float ddy = (float)(bbox->pMax.y - bbox->pMin.y) / outdimy;

		fprintf(file, "x = ");
		for (int i = 0; i < outdimx-1; i++)
			fprintf(file, "%.2f, ", bbox->pMin.x + ddx * i);
		fprintf(file, "%.2f ;\n", bbox->pMin.x + ddx * outdimx);

		fprintf(file, "y = ");
		for (int i = 0; i < outdimy-1; i++)
			fprintf(file, "%.2f, ", bbox->pMin.y + ddy * i);
		fprintf(file, "%.2f ;\n", bbox->pMin.y + ddy * outdimy);
		
		fprintf(file, "time = ");
		for (float cur = 0; cur < time; cur += (float)timestep)
			fprintf(file, "%.2f, ", cur);
		fprintf(file, "%.2f ;\n", time);

		fprintf(file, "u = \n");

		fclose(file);
	}

	static void OutputNetCDF3D_U(const char *outputPath, Vec3D *v, double *T, int dimx, int dimy, int dimz, bool finish)
	{
		FILE *file = NULL;
		fopen_s(&file, outputPath, "a");
		
			for (int i = 0; i < dimx; i++)
			{
				for (int j = 0; j < dimy; j++)
				{
					for (int k = 0; k < dimz; k++)
					{
						fprintf(file, "%.3f", v[i * dimy * dimz + j * dimz + k].x);
						if (finish && (i == dimx-1) && (j == dimy-1) && (k == dimz-1)) fprintf(file, " ; ");
							else fprintf(file, ", ");
					}
					fprintf(file, "\n");
				}
				fprintf(file, "\n");
			}

		if (finish) fprintf(file, "}");
		fclose(file);
	}

	static void OutputNetCDF2D_U(const char *outputPath, Vec2D *v, double *T, int dimx, int dimy, bool finish)
	{
		FILE *file = NULL;
		fopen_s(&file, outputPath, "a");
		
			for (int i = 0; i < dimx; i++)
			{
				for (int j = 0; j < dimy; j++)
				{	
					fprintf(file, "%.3f", v[i * dimy + j].x);
					if (finish && (i == dimx-1) && (j == dimy-1)) fprintf(file, " ; ");
						else fprintf(file, ", ");	
				}
				fprintf(file, "\n");
			}

		if (finish) fprintf(file, "}");
		fclose(file);
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
					v[i * dimy + j].x = (FTYPE)vx;
					v[i * dimy + j].y = (FTYPE)vy;
					T[i * dimy + j] = (FTYPE)t;
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
		if (file) fclose(file);
	}

	static void ReadLine(FILE *file, char* str, int maxlength)
	{
		char c;
		int i = 0;
		memset(str, 0, maxlength);
		for (i=0; i<maxlength-1; i++)
		{
			int s = fread(&c, 1, 1, file);
			if (c!='\n' && s>0)
				str[i] = c;
			else
				break;
		}
	}

	static void ReadPoint2D(FILE *file, Vec2D &p)
	{
		std::string str = "";
		char c;	

		// read line
		fscanf_s(file, "%c", &c);
		while (c == '\n' || c == ' ') fscanf_s(file, "%c", &c);
		while (c != '\n')
		{
			str += c;	
			fscanf_s(file, "%c", &c);
		}

		// replace ',' with '.' if necessary
		std::string::size_type pos = 0;
		std::string::size_type found;
		while ((found = str.find(',', pos)) != std::string::npos)
		{
			str.replace(found, 1, 1, '.');
			pos = found;
		}

		// separate 2 values
		pos = 0;
		found = str.find(' ', pos);
		std::string s1 = str.substr(0, found);
		std::string s2 = str.substr(found+1, std::string::npos);

		// convert to fp
		p.x = (FTYPE)atof(s1.c_str());
		p.y = (FTYPE)atof(s2.c_str());
	}

	static void ReadPoint3D(FILE *file, Vec3D &p)
	{
		std::string str[3];
		str[0] = str[1] = str[2] = "";
		char c;	

		// read line
		fscanf_s(file, "%c", &c);
		while (c == '\n' || c == ' ') fscanf_s(file, "%c", &c);
		
		for (int i = 0; i < 3; i++)
		{
			while (c != '\n' && c != ' ')
			{
				str[i] += c;	
				fscanf_s(file, "%c", &c);
			}

			// replace ',' with '.' if necessary
			std::string::size_type found;
			if ((found = str[i].find(',', 0)) != std::string::npos)
				str[i].replace(found, 1, 1, '.');

			while (c == ' ') fscanf_s(file, "%c", &c);
		}

		// convert to doubles
		p.x = (FTYPE)atof(str[0].c_str());
		p.y = (FTYPE)atof(str[1].c_str());
		p.z = (FTYPE)atof(str[2].c_str());
	}

	static int ExtractInt(char* str)
	{
		int result = 0;
		for(int i=0; i>=0; i++)
		{
			char c = str[i];
			if (c >= '0' && c<= '9')
				result = result*10 + c - '0';
			if (c == 0)
				return result;
		}
		return result;
	}

	static void LoadProject(char *proj, char* inputPath, char* fieldPath, char* outputPath, char* configPath)
	{
		char projectPath[MAX_STR_SIZE];

		char t1[MAX_STR_SIZE];
		char t2[MAX_STR_SIZE];
		char t3[MAX_STR_SIZE];
		char t4[MAX_STR_SIZE];

		FindFile(projectPath, proj);
		FILE* prj = NULL;
		fopen_s(&prj, projectPath, "r");

		ReadLine(prj, t1, MAX_STR_SIZE);
		ReadLine(prj, t2, MAX_STR_SIZE);
		ReadLine(prj, t3, MAX_STR_SIZE);
		ReadLine(prj, t4, MAX_STR_SIZE);

		if (t4[0] != 0)
		{
			FindFile(inputPath, t1);
			FindFile(fieldPath, t2);
			FindFile(outputPath, t3, false);
			FindFile(configPath, t4);
		}
		else
		{
			FindFile(inputPath, t1);
			FindFile(outputPath, t2, false);
			FindFile(configPath, t3);
			sprintf_s(fieldPath, MAX_STR_SIZE, "");
		}

		fclose(prj);
	}

	static void ExtendFileName(char* src, char* dest, char* add)
	{
		int l = strlen(src);
		char temp[MAX_STR_SIZE];

		strcpy_s(temp, src);
		_strrev(temp);
		char* name = strchr(temp, '.') + 1;
		_strrev(name);
		int ln = strlen(name);
		temp[l-ln-1] = 0;
		_strrev(temp);

		sprintf_s(dest, MAX_STR_SIZE, "%s%s.%s", name, add, temp);
	}
}
