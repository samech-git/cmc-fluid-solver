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

#include <string>

#define MAX_STR_SIZE	255

namespace Common
{
	enum solver { Explicit, ADI, Stable, _unknownSolver };
	enum dimension { _2D, _3D, _unknownDim };
	enum inFormat { Shape2D, Shape3D, SeaNetCDF, _unknownInFmt };
	enum outFormat { NetCDF, MultiVox, _unknownOutFmt };
	
	class Config
	{
	public:

		static dimension problem_dim;

		// input format
		static inFormat in_fmt;

		// grid size
		static double dx, dy, dz;

		// fluid parameters
		static bool useNormalizedParams;
		static double viscosity, density;
		static double Re, Pr, lambda;		

		// boundary conditions
		static bool bc_noslip;
		static double bc_strength;			// [0..1] if 0 - noslip, if 1 - slip
		static Vec3D bc_inV;
		static double bc_inT;

		// depth of our 3D object (along Z direction)
		static double depth;				

		// thermodynamic params
		static double R_specific, k, cv, baseT;		 

		// time params
		static int cycles, time_steps, out_time_steps;
		static double frame_time;

		// output grid
		static outFormat out_fmt;
		static int outdimx, outdimy, outdimz;
		static vector<string> out_vars;

		// solver params
		static solver solverID;		
		static int num_global, num_local;

		Config()
		{
			// set default
			R_specific = 461.495;		// water,	287.058	for air (gas constant)
			k = 0.6;					// water (thermal conductivity)
			cv = 4200.0;				// water (specific heat capacity at constant volume)
			baseT = 1.0;				// normalized

			bc_noslip = true;	
			bc_strength = 0.5;
			bc_inV = Vec3D(0.0f, 0.0f, 0.0f);
			bc_inT = baseT;

			useNormalizedParams = false;
			viscosity = 0.05;
			density = 1000.0;
			Re = Pr = lambda = -1;

			cycles = 1;
			time_steps = 50;
			out_time_steps = 10;
			outdimx = outdimy = outdimz = 50;
			out_vars.clear();

			num_global = 2;
			num_local = 1;

			// must specify 
			problem_dim = _unknownDim;
			in_fmt = _unknownInFmt;
			out_fmt = _unknownOutFmt;
			solverID = _unknownSolver;
			frame_time = -1;
			dx = -1;
			dy = -1;
			dz = -1;
			depth = -1;
		}

		static void ReadDouble(FILE *file, double &value)
		{
			float f;
			ReadFloat(file, f);
			value = (double)f;
		}

		static void ReadFloat(FILE *file, float &value)
		{
			float f = 0.0f;
			fscanf_s(file, "%f", &f);
			value = f;
		}

		static void ReadInt(FILE *file, int &value)
		{
			fscanf_s(file, "%i", &value);
		}

		static void ReadSolver(FILE *file)
		{
			char solverStr[MAX_STR_SIZE];
			fscanf_s(file, "%s", solverStr, MAX_STR_SIZE);
			if (!strcmp(solverStr, "Explicit")) solverID = Explicit;
			if (!strcmp(solverStr, "ADI"))		solverID = ADI;
			if (!strcmp(solverStr, "Stable"))	solverID = Stable;
		}

		static void ReadBC(FILE *file)
		{
			char bcStr[MAX_STR_SIZE];
			fscanf_s(file, "%s", bcStr, MAX_STR_SIZE);
			if (!strcmp(bcStr, "NoSlip")) bc_noslip = true;
				else bc_noslip = false;
		}

		static void ReadDim(FILE *file)
		{
			char dimStr[MAX_STR_SIZE];
			fscanf_s(file, "%s", dimStr, MAX_STR_SIZE);
			if (!strcmp(dimStr, "2D")) problem_dim = _2D;
				else problem_dim = _3D;
		}

		static void ReadInFormat(FILE *file)
		{
			char dimStr[MAX_STR_SIZE];
			fscanf_s(file, "%s", dimStr, MAX_STR_SIZE);
			if (!strcmp(dimStr, "Shape2D")) in_fmt = Shape2D;
			else if (!strcmp(dimStr, "Shape3D")) in_fmt = Shape3D;
			else if (!strcmp(dimStr, "SeaNetCDF")) in_fmt = SeaNetCDF;
		}

		static void ReadOutFormat(FILE *file)
		{
			char dimStr[MAX_STR_SIZE];
			fscanf_s(file, "%s", dimStr, MAX_STR_SIZE);
			if (!strcmp(dimStr, "NetCDF")) out_fmt = NetCDF;
				else out_fmt = MultiVox;
		}

		static void ReadVars(FILE *file)
		{
			int num;
			fscanf_s(file, "%i", &num);
			for( int i = 0; i < num; i++ ) {
				char str[MAX_STR_SIZE];
				fscanf_s(file, "%s", str, MAX_STR_SIZE);
				out_vars.push_back( str );
			}
		}

		static void LoadFromFile(char *filename)
		{
			FILE *file = NULL;
			fopen_s(&file, filename, "r");

			if (file == NULL) { printf("cannot open config file!\n"); exit(0); }

			char str[MAX_STR_SIZE];
			while (!feof(file))
			{
				int ret = fscanf_s(file, "%s", str, MAX_STR_SIZE);
				if (ret <= 0) break;

				if (!strcmp(str, "dimension")) ReadDim(file);
				if (!strcmp(str, "in_fmt")) ReadInFormat(file);

				if (!strcmp(str, "viscosity")) ReadDouble(file, viscosity);
				if (!strcmp(str, "density")) ReadDouble(file, density);

				if (!strcmp(str, "Re")) { useNormalizedParams = true; ReadDouble(file, Re); }
				if (!strcmp(str, "Pr")) { useNormalizedParams = true; ReadDouble(file, Pr); }
				if (!strcmp(str, "lambda")) { useNormalizedParams = true; ReadDouble(file, lambda); }

				if (!strcmp(str, "bc_type")) ReadBC(file);
				if (!strcmp(str, "bc_strenght")) ReadDouble(file, bc_strength);
				//if (!strcmp(str, "bc_initv")) { ReadDouble(file, bc_inV.x); ReadDouble(file, bc_inV.y); ReadDouble(file, bc_inV.z); }
				if (!strcmp(str, "bc_initv")) { ReadFloat(file, bc_inV.x); ReadFloat(file, bc_inV.y); ReadFloat(file, bc_inV.z); }
				if (!strcmp(str, "bc_initT")) { ReadDouble(file, bc_inT); }

				if (!strcmp(str, "grid_dx")) ReadDouble(file, dx);
				if (!strcmp(str, "grid_dy")) ReadDouble(file, dy);
				if (!strcmp(str, "grid_dz")) ReadDouble(file, dz);

				if (!strcmp(str, "cycles")) ReadInt(file, cycles);
				if (!strcmp(str, "frame_time")) ReadDouble(file, frame_time);
				if (!strcmp(str, "time_steps")) ReadInt(file, time_steps);

				if (!strcmp(str, "out_vars")) ReadVars(file);
				if (!strcmp(str, "out_time_steps")) ReadInt(file, out_time_steps);
				if (!strcmp(str, "out_gridx")) ReadInt(file, outdimx);
				if (!strcmp(str, "out_gridy")) ReadInt(file, outdimy);
				if (!strcmp(str, "out_gridz")) ReadInt(file, outdimz);
				if (!strcmp(str, "out_fmt")) ReadOutFormat(file);
				
				if (!strcmp(str, "depth")) ReadDouble(file, depth);		

				if (!strcmp(str, "solver")) ReadSolver(file);
				if (!strcmp(str, "num_global")) ReadInt(file, num_global);
				if (!strcmp(str, "num_local")) ReadInt(file, num_local);
			}	

			fclose(file);

			// checking		
			if (problem_dim == _unknownDim) { printf("must specify problem dimension!\n"); exit(0); }
			if (solverID == _unknownSolver) { printf("must specify solver!\n"); exit(0); }
			if (out_fmt == _unknownOutFmt) { printf("must specify output format!\n"); exit(0); }
			
			if (frame_time < 0 && in_fmt == SeaNetCDF) { printf("must specify frame time!\n"); exit(0); }
			if (dx < 0) { printf("cannot find dx!\n"); exit(0); }
			if (dy < 0) { printf("cannot find dy!\n"); exit(0); }
			
			if (problem_dim == _2D) in_fmt = Shape2D;
			if (problem_dim == _3D)
			{
				if (out_vars.empty()) { printf("must output at least 1 var!\n"); exit(0); }
				if (in_fmt == _unknownInFmt) { printf("must specify input format!\n"); exit(0); }
				if (dz < 0) { printf("cannot find dz!\n"); exit(0); }
				if (in_fmt == Shape2D)
				{
					if (depth < 0) { printf("cannot find depth!\n"); exit(0); }
				}
				if (out_fmt == MultiVox) { printf("MultiVox output format is not supported for 3D modes\n"); exit(0); }
			}
			if (useNormalizedParams && (Re < 0 || Pr < 0 || lambda < 0)) { printf("must specify Re, Pr and lambda!\n"); exit(0); }
		}
	};

	dimension Config::problem_dim;
	inFormat Config::in_fmt;
	outFormat Config::out_fmt;

	double Config::dx, Config::dy, Config::dz;

	bool Config::useNormalizedParams;
	double Config::viscosity, Config::density;
	double Config::Re, Config::Pr, Config::lambda;	

	bool Config::bc_noslip;
	double Config::bc_strength;
	Vec3D Config::bc_inV;
	double Config::bc_inT;

	double Config::depth;

	double Config::R_specific, Config::k, Config::cv, Config::baseT;		 

	int Config::cycles, Config::time_steps, Config::out_time_steps;
	double Config::frame_time;

	int Config::outdimx, Config::outdimy, Config::outdimz;
	vector<string> Config::out_vars;

	solver Config::solverID;		
	int Config::num_global, Config::num_local;
}
