#pragma once

#include "../Common/Common.h"

#include <stdlib.h>
#include <math.h>

#include <vector>
#include <string>

#define BBOX_OFFSET			0.03
#define GRID_SCALE_FACTOR	0.001f
#define MAX_STR_SIZE		255

using namespace std;
using namespace Common;

namespace FluidSolver2D
{
	enum CellType { IN, OUT, BOUND, VALVE };
	enum CondType { NONE, NOSLIP, FREE };

	struct CondData2D
	{
		CellType cell;
		CondType type;
		Vec2D vel;
		double T;

		CondData2D() : type(NONE), vel(Vec2D(0.0, 0.0)), T(0.0) { }
		CondData2D(CondType _type, CellType _cell, Vec2D _vel, double _T) : type(_type), vel(_vel), T(_T), cell(_cell) { }
	};

	struct Grid2D
	{
		Grid2D(double _dx, double _dy, double _startT, bool _bc_slip, double _bc_strength);
		Grid2D(Grid2D &grid);
		~Grid2D();

		int dimx, dimy;
		double dx, dy;

		CellType GetType(int x, int y);
		CondData2D GetData(int x, int y);
		void SetFieldData(int x, int y, CondData2D d);
		
		void Prepare(int frame, double substep);
		void Prepare(double time);
		double GetCycleLenght();
		int GetFramesNum();
		int GetFrame(double time);
		float GetLayerTime(double t);

		bool LoadFromFile(char *filename, char *fieldname);
		void TestPrint(char *filename);

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

		void Init();

		VecTN GetTangentNormal(Vec2D vector, Vec2D orientation);
		Vec2D GetBounfVelocity(int x, int y);
		void RasterLine(Point2D p1, Point2D p2, Vec2D v1, Vec2D v2, CellType color);
		void RasterField(Field2D field);
		void FloodFill(CellType color);
		
		void Build(FrameInfo2D frame);
		
		void SetType(int x, int y, CellType t);
		void SetData(int x, int y, CondData2D d);
	};
}