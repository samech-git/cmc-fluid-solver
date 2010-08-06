#pragma once

#include <string>

#define MAX_STR_SIZE	255

namespace Common
{
	enum solver { Explicit, ADI, Stable };
	enum dimension { _2D, _3D, unknown };
	
	class Config
	{
	public:

		static dimension problem_dim;

		// grid size
		static double dx, dy, dz;

		// fluid parameters
		static bool useNormalizedParams;
		static double viscosity, density;
		static double Re, Pr, lambda;		

		// boundary conditions
		static bool bc_noslip;
		static double bc_strength;			// [0..1] if 0 - noslip, if 1 - slip

		// depth of our 3D object (along Z direction)
		static double depth;				

		// thermodynamic params
		static double R_specific, k, cv, startT;		 

		// animation params
		static int cycles, calc_subframes, out_subframes;

		// output grid
		static int outdimx, outdimy, outdimz;

		// solver params
		static solver solverID;		
		static int num_global, num_local;

		Config()
		{
			// set default
			R_specific = 461.495;		// water,	287.058	for air (gas constant)
			k = 0.6;					// water (thermal conductivity)
			cv = 4200.0;				// water (specific heat capacity at constant volume)
			startT = 1.0;				// normalized

			bc_noslip = true;	
			bc_strength = 0.5;

			useNormalizedParams = false;
			viscosity = 0.05;
			density = 1000.0;
			Re = Pr = lambda = -1;

			cycles = 1;
			calc_subframes = 50;
			out_subframes = 10;
			outdimx = outdimy = outdimz = 50;

			solverID = Stable;
			num_global = 2;
			num_local = 1;

			// must specify 
			problem_dim = unknown;
			dx = -1;
			dy = -1;
			dz = -1;
			depth = -1;
		}

		static void ReadDouble(FILE *file, double &value)
		{
			float f = 0.0f;
			fscanf_s(file, "%f", &f);
			value = (double)f;
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

				if (!strcmp(str, "viscosity")) ReadDouble(file, viscosity);
				if (!strcmp(str, "density")) ReadDouble(file, density);

				if (!strcmp(str, "Re")) { useNormalizedParams = true; ReadDouble(file, Re); }
				if (!strcmp(str, "Pr")) { useNormalizedParams = true; ReadDouble(file, Pr); }
				if (!strcmp(str, "lambda")) { useNormalizedParams = true; ReadDouble(file, lambda); }

				if (!strcmp(str, "bc_type")) ReadBC(file);
				if (!strcmp(str, "bc_strenght")) ReadDouble(file, bc_strength);

				if (!strcmp(str, "grid_dx")) ReadDouble(file, dx);
				if (!strcmp(str, "grid_dy")) ReadDouble(file, dy);
				if (!strcmp(str, "grid_dz")) ReadDouble(file, dz);

				if (!strcmp(str, "cycles")) ReadInt(file, cycles);
				if (!strcmp(str, "calc_subframes")) ReadInt(file, calc_subframes);

				if (!strcmp(str, "out_subframes")) ReadInt(file, out_subframes);
				if (!strcmp(str, "out_gridx")) ReadInt(file, outdimx);
				if (!strcmp(str, "out_gridy")) ReadInt(file, outdimy);
				if (!strcmp(str, "out_gridz")) ReadInt(file, outdimz);

				if (!strcmp(str, "depth")) ReadDouble(file, depth);		

				if (!strcmp(str, "solver")) ReadSolver(file);
				if (!strcmp(str, "num_global")) ReadInt(file, num_global);
				if (!strcmp(str, "num_local")) ReadInt(file, num_local);
			}	

			fclose(file);

			// checking		
			if (problem_dim == unknown) { printf("must specify problem dimension!\n"); exit(0); }
			if (dx < 0) { printf("cannot find dx!\n"); exit(0); }
			if (dy < 0) { printf("cannot find dy!\n"); exit(0); }
			if (problem_dim == _3D)
			{
				if (dz < 0) { printf("cannot find dz!\n"); exit(0); }
				if (depth < 0) { printf("cannot find depth!\n"); exit(0); }
			}
			if (useNormalizedParams && (Re < 0 || Pr < 0 || lambda < 0)) { printf("must specify Re, Pr and lambda!\n"); exit(0); }
		}
	};

	dimension Config::problem_dim;

	double Config::dx, Config::dy, Config::dz;

	bool Config::useNormalizedParams;
	double Config::viscosity, Config::density;
	double Config::Re, Config::Pr, Config::lambda;	

	bool Config::bc_noslip;

	double Config::bc_strength;

	double Config::depth;

	double Config::R_specific, Config::k, Config::cv, Config::startT;		 

	int Config::cycles, Config::calc_subframes, Config::out_subframes;

	int Config::outdimx, Config::outdimy, Config::outdimz;

	solver Config::solverID;		
	int Config::num_global, Config::num_local;
}