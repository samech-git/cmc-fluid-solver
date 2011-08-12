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
#include "../Common/test_util.h"
//#include <conio.h> // for getch

#ifdef linux
#include <typeinfo>
#endif

using namespace FluidSolver3D;
using namespace Common;

void parse_cmd_params(int argc, char **argv, BackendType &backend, bool &csv, bool &transpose, bool &decompose, bool &align, int &nGPU)
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
		int nGPU = 0;
		parse_cmd_params(argc, argv, backend, csv, transpose, decompose, align, nGPU);

#define SIMPLE_GRID 0
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
#if !SIMPLE_GRID
		Grid3D *grid = NULL;
		if( Config::in_fmt == Shape3D ) 
		{
			grid = new Grid3D(Config::dx, Config::dy, Config::dz, Config::baseT, backend);
			if (pplan->rank() == 0)
				printf("Geometry: 3D polygons\n", grid->dimx, grid->dimy, grid->dimz);
		}
		else if( Config::in_fmt == Shape2D )
		{
			grid = new Grid3D(Config::dx, Config::dy, Config::dz, Config::depth, Config::baseT, backend);
			if (pplan->rank() == 0)
				printf("Geometry: extruded 2D shape\n");
		}
		else
		{
			grid = new Grid3D(Config::dx, Config::dy, Config::dz, Config::baseT, backend, true);
			if (pplan->rank() == 0)
				printf("Geometry: depths from NetCDF\n");
		}

		grid->SetFrameTime( Config::frame_time );
		grid->SetBoundParams( Config::bc_inV, Config::bc_inT );

		printf("Grid options:\n  align %s\n", align ? "ON" : "OFF");
		if (grid->LoadFromFile(inputPath, align))
			if (pplan->rank() == 0)
				printf("Grid = %i x %i x %i\n", grid->dimx, grid->dimy, grid->dimz);
		grid->Prepare(0.0);

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
		//grid->printTypes();
#endif

		//--------------------------------------- For Testing ---------------------------------------
#if SIMPLE_GRID
		Grid3D *grid = NULL;
		grid = new Grid3D(Config::dx, Config::dy, Config::dz, Config::baseT, backend);
		grid->SetFrameTime( Config::frame_time );

		grid->SetFrameTime( Config::frame_time );
		grid->SetBoundParams( Config::bc_inV, Config::bc_inT );
		grid->genRandom();
		printf("Geometry: rectangular shape with random T data:  %d	%d	%d\n", grid->dimx, grid->dimy, grid->dimz);

		grid->Prepare2();

		//grid->printTypes();

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
					printf("Solver options:\n  transpose %s\n  decompose %s\n", transpose ? "ON" : "OFF", decompose ? "ON" : "OFF");
					dynamic_cast<AdiSolver3D*>(solver)->SetOptionsGPU(transpose, decompose);
				}
				break;
		}

		// allocate result arrays			
		int outsize = Config::outdimy * Config::outdimz;
		int noutdimx, outoffset;
		pplan->getEven1D(noutdimx, outoffset, Config::outdimx);
		outsize *= (pplan->rank()==0)? Config::outdimx : noutdimx;
		Vec3D *resVel = new Vec3D[outsize];
		double *resT = new double[outsize];

		solver->Init(backend, csv, grid, *params);
		
		solver->UpdateBoundaries();
		solver->TimeStep((FTYPE)0.1, Config::num_global, Config::num_local);
		solver->SetGridBoundaries();
#endif
		//-------------------------------------------------------------------------------------------

#if !SIMPLE_GRID
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
						printf("Solver options:\n  transpose %s\n  decompose %s\n", transpose ? "ON" : "OFF", decompose ? "ON" : "OFF");
					dynamic_cast<AdiSolver3D*>(solver)->SetOptionsGPU(transpose, decompose);
				}
				break;
		}
		solver->Init(backend, csv, grid, *params);
		int startFrame = 0;

		int frames = grid->GetFramesNum();
		double length = grid->GetCycleLength();
		double dt = length / (frames * Config::time_steps);
		double finaltime = length * Config::cycles;
		if (pplan->rank() == 0)
		{
			// create file and output header
			sprintf_s(outputPath, MAX_STR_SIZE, "%s_res.nc", argv[2]);
			BBox3D *bbox = NULL;
			if( Config::in_fmt == Shape2D ) bbox = new BBox3D( grid->GetGrid2D()->bbox, (float)Config::depth );
				else bbox = &grid->GetBBox();
			OutputNetCDF3D_header(outputPath, bbox, grid->GetDepthInfo(), dt * Config::out_time_steps, finaltime, Config::outdimx, Config::outdimy, Config::outdimz, Config::out_vars, Config::in_fmt == SeaNetCDF );
			fflush(stdout);
		}

			// allocate result arrays			
			int outsize = Config::outdimy * Config::outdimz;
			int noutdimx, outoffset;
			pplan->getEven1D(noutdimx, outoffset, Config::outdimx);
			outsize *= (pplan->rank()==0)? Config::outdimx : noutdimx;
			Vec3D *resVel = new Vec3D[outsize];
			double *resT = new double[outsize];
//
//		//------------------------------------------ Solving ------------------------------------------
		cpu_timer timer;
		timer.start();
		int lastframe = -1;
		int out_layer = 0;
		double t = dt;
		for (int i=0; t < finaltime; t+=dt, i++)
		{
			int currentframe = grid->GetFrame(t);
			float layer_time = grid->GetLayerTime(t);

			if (currentframe != lastframe)
			{
				lastframe = currentframe;
				i = 0;
			}

			grid->Prepare(t);

			solver->UpdateBoundaries();
			solver->TimeStep((FTYPE)dt, Config::num_global, Config::num_local);
			solver->SetGridBoundaries();

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
		if (pplan->rank() == 0)
		{
			printf("\nTotal time: %.2f sec\n", timer.elapsed_sec());
			OutputSliceResult("last_zslice.txt", 18, resVel, resT, Config::outdimx, Config::outdimy, Config::outdimz, (float)dt * Config::out_time_steps);
		}
#endif

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
		//getch();
#ifdef __PARA
		MPI_Abort(MPI_COMM_WORLD, -1);
#endif
		return -1;
	}
	fflush(stdout);
	//getch();
	return 0;
}
