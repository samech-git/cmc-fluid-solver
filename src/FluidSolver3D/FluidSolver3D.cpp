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

#include "FluidSolver3D.h"

#ifdef linux
#include <typeinfo>
#endif

using namespace FluidSolver3D;
using namespace Common;

void parse_cmd_params(int argc, char **argv, BackendType &backend, bool &csv, bool &transpose, bool &decompose, bool &align, int &nGPU, bool &blocking, int &nBlockZ)
{
	for( int i = 4; i < argc; i++ )
	{
		if( !strcmp(argv[i], "GPU") )
		{
			backend = GPU;
			if ( i == argc-1)
				nGPU = 0;
			else
				nGPU = atoi(argv[i+1]);
		}
		if( !strcmp(argv[i], "blocking") )
			{
				blocking = true;
				if ( i == argc-1)
					nBlockZ = 1;
				else
					nBlockZ = atoi(argv[i+1]);
		}
		if( !strcmp(argv[i], "CSV") ) csv = true;
		if( !strcmp(argv[i], "transpose") ) transpose = true;
		if( !strcmp(argv[i], "decompose") ) decompose = true;
		if( !strcmp(argv[i], "align") ) align = true;
	}
}

int main(int argc, char **argv)
{
	try
	{
#ifdef __PARA
		MPI_Init(&argc, &argv);
#endif
		BackendType backend = CPU;
		PARAplan *pplan = PARAplan::Instance();		
		bool csv = false;
		bool transpose = false;
		bool decompose = false;
		bool align = false;
		bool useBlocks = false;
		int nBlockZ = 1;
		int nGPU = 0;
		parse_cmd_params(argc, argv, backend, csv, transpose, decompose, align, nGPU, useBlocks, nBlockZ);

		pplan->init(backend);
		if( backend == CPU )
		{
			if (pplan->size() > 1)
				throw runtime_error("MPI is not supported for backend = 'CPU'");
#ifdef _OPENMP
			printf("Using OpenMP: num_proc = %i\n", omp_get_num_procs());
#else
			printf("Using single CPU\n");
#endif
		}
		else
		{
			if (nGPU != 0)
				pplan->setGPUnum(nGPU);

			cudaDeviceProp deviceProp;
#if MGPU_EMU
			printf("Multi GPU emulation is ON\n");
#endif
			for (int i = 0; i < pplan->gpuNum(); i++)
			{
#if MGPU_EMU
				cudaGetDeviceProperties(&deviceProp, DEFAULT_DEVICE);
				printf("Node(%d): Using GPU %d: %s\n", pplan->rank(), DEFAULT_DEVICE, deviceProp.name);
#else
				cudaGetDeviceProperties(&deviceProp, i);
				printf("Node(%d): Using GPU %d: %s\n", pplan->rank(), i, deviceProp.name);
				fflush(stdout);
#endif
			}
		}

		if (pplan->rank() == 0)
			printf("%s precision computations\n", (typeid(FTYPE) == typeid(float)) ?  "Single" : "Double");

		char inputPath[MAX_STR_SIZE];
		char configPath[MAX_STR_SIZE];
		char outputPath[MAX_STR_SIZE];
		char gridPath[MAX_STR_SIZE];

		FindFile(inputPath, argv[1]);
		FindFile(configPath, argv[3]);
	
		Config::Config();
		Config::LoadFromFile(configPath);
		//--------------------------------------- Initializing ---------------------------------------
		Grid3D *grid = NULL;
		SplitType split_type = EVEN_X; //EVEN_X, EVEN_SEGMENTS or EVEN_VOLUME

		if( Config::in_fmt == Shape3D ) 
		{
			grid = new Grid3D(Config::dx, Config::dy, Config::dz, Config::baseT, backend, false, split_type);
			if (pplan->rank() == 0)
				printf("Geometry: 3D polygons\n", grid->dimx, grid->dimy, grid->dimz);
		}
		else if( Config::in_fmt == Shape2D )
		{
			grid = new Grid3D(Config::dx, Config::dy, Config::dz, Config::depth, Config::baseT, backend, false, split_type);
			if (pplan->rank() == 0)
				printf("Geometry: extruded 2D shape\n");
		}
		else
		{
			grid = new Grid3D(Config::dx, Config::dy, Config::dz, Config::baseT, backend, true, split_type);
			if (pplan->rank() == 0)
				printf("Geometry: depths from NetCDF\n");
		}

		grid->SetFrameTime( Config::frame_time );
		grid->SetBoundParams( Config::bc_inV, Config::bc_inT );

		printf("Grid options:\n  align %s\n", align ? "ON" : "OFF");
		if (grid->LoadFromFile(inputPath, align))
			if (pplan->rank() == 0)
				printf("Grid = %i x %i x %i\n", grid->dimx, grid->dimy, grid->dimz);
		grid->Prepare_CPU(0.0);
		grid->Split();
		grid->Init_GPU();
		if (pplan->rank() == 0)
		{
			sprintf_s(gridPath, "%s_grid_3d", argv[2]);
			grid->OutputImage(gridPath);

			if( grid->GetGrid2D() != NULL )
			{
				sprintf_s(gridPath, "%s_grid_2d.bmp", argv[2]);
				grid->GetGrid2D()->OutputImage(gridPath);
			}	
		}

		//calculate volume:
		double indsidePoints = 0.;
		for (int i = 0; i < grid->dimx; i++)
			for (int j = 0; j < grid->dimy; j++)
				for (int k = 0; k < grid->dimz; k++)
					if (grid->GetType(i,j,k) == NODE_IN)
						indsidePoints += 1.0;
		if (pplan->rank() == 0)
			printf("NODE_IN points = %f of total %f, volume = %f\n", indsidePoints, double(grid->dimx) * grid->dimy * grid->dimz, indsidePoints * grid->dx * grid->dy * grid->dz);

		FluidParams *params;
		if (Config::useNormalizedParams) params = new FluidParams(Config::Re, Config::Pr, Config::lambda);
			else params = new FluidParams(Config::viscosity, Config::density, Config::R_specific, Config::k, Config::cv);

		Solver3D *solver;
		switch (Config::solverID)
		{
			case Explicit: printf("Explicit solver is not implemented yet!\n"); break;
			case Stable: printf("Stable solver is not implemented yet!\n"); break;
			case ADI: 
				solver = new AdiSolver3D(); 
				if( backend == GPU ) 
				{
					if (pplan->rank() == 0)
						printf("Solver options:\n  transpose %s\n  decompose %s\n  number of blocks %d\n", transpose ? "ON" : "OFF", decompose ? "ON" : "OFF", nBlockZ);
					dynamic_cast<AdiSolver3D*>(solver)->SetOptionsGPU(transpose, decompose);
				}
				break;
		}
		solver->Init(backend, csv, grid, *params, useBlocks, nBlockZ);

		int startFrame = 0;

		int frames = grid->GetFramesNum();
		double length = grid->GetCycleLength();
		double dt = length / (frames * Config::time_steps);
		double finaltime = length * Config::cycles;
		if (pplan->rank() == 0)
		{
			// create file and output header
			BBox3D *bbox = NULL;
			sprintf_s(outputPath, MAX_STR_SIZE, "%s_res.nc", argv[2]);
			if( Config::in_fmt == Shape2D ) bbox = new BBox3D( grid->GetGrid2D()->bbox, (float)Config::depth );
				else bbox = &grid->GetBBox();
			OutputNetCDF3D_header(outputPath, bbox, grid->GetDepthInfo(), dt * Config::out_time_steps, finaltime, Config::outdimx, Config::outdimy, Config::outdimz, Config::out_vars, Config::in_fmt == SeaNetCDF );
			fflush(stdout);
		}

			// allocate result arrays			
			int outsize = Config::outdimy * Config::outdimz;
			int noutdimx, outoffset;
			pplan->get1D(noutdimx, outoffset, Config::outdimx);
			outsize *= (pplan->rank()==0)? Config::outdimx : noutdimx;
			Vec3D *resVel = new Vec3D[outsize];
			double *resT = new double[outsize];

		//------------------------------------------ Solving ------------------------------------------
		cpu_timer timer;
		timer.start();
		int lastframe = -1;
		int out_layer = 0;
		double t = dt;
		dynamic_cast<AdiSolver3D*>(solver)->CreateSegments();
		grid->Prepare(0);
		for (int i=0; t < finaltime; t+=dt, i++)
		{
			int currentframe = grid->GetFrame(t);
			float layer_time = grid->GetLayerTime(t);

			if (currentframe != lastframe)
			{
				lastframe = currentframe;
				i = 0;
			}

			//grid->Prepare(t);

			/*if (i == 0)
				solver->debug(true);*/			
			solver->UpdateBoundaries(); // needs this since cur gets overwritten (do not call CreateSegments, so it is quite cheap)
			solver->TimeStep((FTYPE)dt, Config::num_global, Config::num_local, (i%10 == 0) || (t + dt >= finaltime));			
			// solver->SetGridBoundaries();
			//if (i == 0)
			//	solver->debug(false);

			timer.stop();

			PrintTimeStepInfo(currentframe, i, t, finaltime, timer.elapsed_sec());

			if ((i % Config::out_time_steps) == 0)
			{
				float dur = (float)dt * Config::out_time_steps;
				if (dur > layer_time) dur = layer_time;
				solver->GetLayer(resVel, resT, Config::outdimx, Config::outdimy, Config::outdimz);
				if (pplan->rank() == 0)
				{
					OutputNetCDF3D_layer(outputPath, resVel, resT, out_layer, Config::outdimx, Config::outdimy, Config::outdimz, Config::out_vars);
				}
				out_layer++;
			}
		}
		timer.stop();

		delete solver;
		delete [] resVel;
		delete [] resT;

		delete grid;
		delete pplan;
	}
	catch (std::exception& e)
	{
		fprintf(stderr, "\n\nCaught exception:\n");
		fprintf(stderr, "%s\n", e.what());
		fprintf (stderr, "\nTerminating...\n");
		fflush(stdout);
		fflush(stderr);
#ifdef __PARA
		MPI_Abort(MPI_COMM_WORLD, -1);
#endif
		return -1;
	}
	fflush(stdout);
	return 0;
}
