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

#ifdef _WIN32
#include "..\Common\Geometry.h"
#include "..\Common\IO.h"
#include "..\Common\GPUplan.h"
#include "..\Common\PARAplan.h"

#include "..\FluidSolver2D\Grid2D.h"

#elif __unix__
#include "../Common/Geometry.h"
#include "../Common/IO.h"
#include "../Common/GPUplan.h"
#include "../Common/PARAplan.h"

#include "../FluidSolver2D/Grid2D.h"
#endif

#include <cuda_runtime.h>
#include <netcdf.h>

#define TRANSPOSE_OPT 1

using namespace Common;

#define MAX_SEGS_PER_ROW	2

namespace FluidSolver3D
{

	enum SegmentType // Segments get split up along X direction in multiGPU code
	{ 
		BOUND, 
		BOUND_START, 
		BOUND_END,
		UNBOUND
	};

	enum SplitType // Segments get split up along X direction in multiGPU code
	{ 
		EVEN_X, 
		EVEN_SEGMENTS,
		EVEN_VOLUME
	};

	struct Segment3D
	{
		int posx, posy, posz;
		int endx, endy, endz;
		int size;
		bool skipX; // segment alongX to be skipped 
		DirType dir; 
		SegmentType type;
	};

	struct Node
	{
		NodeType type;
		BCtype bc_vel, bc_temp;
		Vec3D v;
		FTYPE T;

		void SetBound(BCtype _bc_vel, BCtype _bc_temp, Vec3D _v, FTYPE _T)
		{
			type = NODE_BOUND;
			bc_vel = _bc_vel;
			bc_temp = _bc_temp;
			v = _v;
			T = _T;
		}
	};

	struct NodesBoundary3D
	{
		Node first, last; // first and last node of a 3D segment
	};

	struct Grid3D
	{
		int dimx, dimy, dimz;
		double dx, dy, dz;

		double baseT;
		Vec3D bcInVel;
		double bcInT;

		Grid3D(double _dx, double _dy, double _dz, double _depth, double _baseT, BackendType _backend, bool useNetCDF = false, SplitType = EVEN_X);			// 2D shape with constant depth
		Grid3D(double _dx, double _dy, double _dz, double _baseT, BackendType _backend, bool useNetCDF = false, SplitType = EVEN_X);						// 3D shape, polygons
		~Grid3D();

		BBox3D GetBBox();

		void GenerateListSegments(int &numSeg, Segment3D *h_list, int dim1, int dim2, int dim3, DirType dir, int nblockZ);
		void GenerateGridBoundaries(NodesBoundary3D *node_list, int numSeg, Segment3D *h_list, bool transposed);
		void SplitSegments_X(int *splitting);

		NodeType GetType(int i, int j, int k);
		BCtype GetBC_vel(int i, int j, int k);
		BCtype GetBC_temp(int i, int j, int k);
		Vec3D GetVel(int i, int j, int k);
		FTYPE GetT(int i, int j, int k);

		void SetType(int i, int j, int k, NodeType type);
		void SetData(int i, int j, int k, BCtype bc_vel, BCtype bc_T, const Vec3D &vel, FTYPE T);

		// return all nodes info as an array
		Node *GetNodesCPU();							 
		NodeType **GetTypesGPU(bool transposed = false);

		void SetNodeVel(int i, int j, int k, Vec3D new_v);

		// frame stuff
		double GetFrameTime();
		int GetFramesNum();
		double GetCycleLength();
		int GetFrame(double time);
		float GetLayerTime(double time);
		void SetFrameTime(double time);
		void SetBoundParams(const Vec3D &vec, const double &temp);

		FluidSolver2D::Grid2D *GetGrid2D();
		DepthInfo3D *GetDepthInfo();
		
		bool LoadFromFile(char *filename, bool align = false);		
		void Prepare(double time);		
		void Prepare_CPU(double time);
		void Init_GPU();
		void Split();

		void TestPrint(char *filename);
		void OutputImage(char *filename_base);

	protected:
		BackendType backend;

		Node*		nodes;		// all grid nodes

		NodeType** d_types; // node types stored on multiple GPUs
		NodeType** d_typesT; // transposed node types on multiple GPU

		bool use3Dshape;
		bool useNetCDF;

		BBox3D bbox;

		// input data
		FrameInfo3D* frames;
		int num_frames;
		double frame_time;
		SplitType split_type;

		// static by now
		DepthInfo3D *depthInfo;

		// support for different input formats
		bool LoadNetCDF(char *filename, bool align);
		bool Load3DShape(char *filename, bool align);

		// helper functions for 3D shape update
		void Init(bool align);
		void Prepare_GPU();
		void Prepare3D_Shape(double time);
		void Prepare3D_NetCDF(double time);

		void ComputeSubframeInfo(int frame, FTYPE substep, FrameInfo3D &res);
		void Build(FrameInfo3D &frame);
		
		void RasterPolygon(Vec3D p1, Vec3D p2, Vec3D p3, Vec3D v1, Vec3D v2, Vec3D v3, NodeType color);
		void RasterLine(Vec3D p1, Vec3D p2, Vec3D v1, Vec3D v2, NodeType color);
		
		double GetBarycentric(Vec2D v1, Vec2D v2, Vec2D v0, Vec2D p);
		void ProjectPointOnPolygon(DirType dir, int i, int j, Vec2D testp, Vec3D n, FTYPE d);
		
		void FloodFill(int start[3], NodeType color, int n, const int *neighborPos);
		
		// helper function for 2D shape update
		void Prepare2D(double time);

		FluidSolver2D::Grid2D *grid2D;		// 2D helper grid for borders
		double depth;						// depth
		int active_dimz;					// dimz before align
	};
}
