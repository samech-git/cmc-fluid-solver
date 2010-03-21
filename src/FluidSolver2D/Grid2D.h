#pragma once

#include "Common.h"

#include <stdlib.h>
#include <math.h>

#include <vector>
#include <string>

#define MAX_STR_SIZE	255
#define BBOX_OFFSET		0.03

using namespace std;

namespace FluidSolver
{
	enum CellType { IN, OUT, BOUND, VALVE };

	struct Shape
	{
		Point2D* Points;
		Vec2D* Velocities;
		int NumPoints;
		bool Active;

		void Init(int num)
		{
			NumPoints = num;
			Points = new Point2D[num];
			Velocities = new Vec2D[num];
		}

		void Dispose()
		{
			delete[] Points;
			delete[] Velocities;
		}
	};

	struct FrameInfo
	{
		Shape* Shapes;
		int NumShapes;
		double Duration;

		void Init(int num)
		{
			NumShapes = num;
			Shapes = new Shape[num];
		}
		void Dispose()
		{
			for (int i=0; i<NumShapes; i++)
				Shapes[i].Dispose();
			delete[] Shapes;
		}
	};

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
		Grid2D(double _dx, double _dy, double _startT, bool _bc_slip);
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

		int LoadFromFile(char *filename);
		void TestPrint();

		BBox2D bbox;		// bounding box
		double startT;		// initial temperature

	protected:
		//---------------------------- Borders ---------------------------------
		FrameInfo* frames;
		int num_frames;

		bool bc_noslip;
		
		void ComputeBorderVelocities(int frame);
		FrameInfo ComputeSubframe(int frame, double substep);
		//----------------------------------------------------------------------

		CondData2D *curData;
		CondData2D *nextData;

		void Init();

		VecTN GetTangentNormal(Vec2D vector, Vec2D orientation);
		Vec2D GetBounfVelocity(int x, int y);
		void RasterLine(Point2D p1, Point2D p2, Vec2D v1, Vec2D v2, CellType color);
		void FloodFill(CellType color);
		
		void BuildBBox(int num_frames, FrameInfo* frames);
		
		void Build(FrameInfo frame);
		
		void SetType(int x, int y, CellType t);
		void SetData(int x, int y, CondData2D d);
	};
}