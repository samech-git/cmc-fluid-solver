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
#elif __unix__
#include "../Common/Geometry.h"
#include "../Common/IO.h"
#endif

#include <stdlib.h>
#include <vector>

#define BBOX_OFFSET			0.03
#define GRID_SCALE_FACTOR	0.001f
#define MAX_STR_SIZE		255

using namespace std;
using namespace Common;

namespace FluidSolver2D
{
	struct CondData2D
	{
		NodeType cell;
		BCtype type;
		
		Vec2D vel;
		FTYPE T;

		CondData2D() : vel(Vec2D(0.0, 0.0)), T(0.0) { }
		CondData2D(BCtype _type, NodeType _cell, Vec2D _vel, FTYPE _T) : type(_type), vel(_vel), T(_T), cell(_cell) { }
	};

	struct Grid2D
	{
		Grid2D(double _dx, double _dy, double _startT, bool _bc_slip, double _bc_strength);
		Grid2D(Grid2D &grid);
		~Grid2D();

		int dimx, dimy;
		double dx, dy;

		inline NodeType GetType(int x, int y);
		CondData2D GetData(int x, int y);
		void SetFieldData(int x, int y, CondData2D d);
		
		void Prepare(int frame, double substep);
		void Prepare(double time);
		double GetCycleLenght();
		int GetFramesNum();
		int GetFrame(double time);
		float GetLayerTime(double t);

		bool LoadFromFile(char *filename, char *fieldname, bool align = false);
		void OutputText(char *filename);
		void OutputImage(char *filename);

		BBox2D bbox;		// bounding box
		double startT;		// initial temperature

	protected:
		//---------------------------- Borders ---------------------------------
		FrameInfo2D* frames;
		int num_frames;

		bool bc_noslip;
		double bc_strength;
		
		void ComputeBorderVelocities(int frame);
		FrameInfo2D ComputeSubframe(int frame, double substep);
		//----------------------------------------------------------------------

		CondData2D *curData;
		CondData2D *nextData;

		void Init(bool align);

		VecTN GetTangentNormal(Vec2D vector, Vec2D orientation);
		Vec2D GetBounfVelocity(int x, int y);
		void RasterLine(Vec2D p1, Vec2D p2, Vec2D v1, Vec2D v2, NodeType color);
		void RasterField(Field2D field);
		void FloodFill(NodeType color);
		
		void Build(FrameInfo2D frame);
		
		void SetType(int x, int y, NodeType t);
		void SetData(int x, int y, CondData2D d);
	};
}
