#pragma once

#include "Common.h"

#include <stdlib.h>
#include <math.h>

#include <vector>
#include <string>

#define INF				1e10
#define MAX_STR_SIZE	255
#define BBOX_OFFSET		3.0

using namespace std;

namespace FluidSolver
{
	enum CellType { IN, OUT, BOUND, VALVE };

	struct Point2D 
	{ 
		double x, y; 

		Point2D() : x(0.0), y(0.0) { }
		Point2D(double _x, double _y) : x(_x), y(_y) { }
	};

	enum CondType { NONE, NOSLIP, FREE };

	struct CondData2D
	{
		CondType type;
		Vec2D vel;
		double T;

		CondData2D() : type(NONE), vel(Vec2D(0.0, 0.0)), T(0.0) { }
		CondData2D(CondType _type, Vec2D _vel, double _T) : type(_type), vel(_vel), T(_T) { }
	};

	struct BBox2D
	{
		Point2D pMin, pMax;

		BBox2D() { Clear(); }
		
		void AddPoint(Point2D p)
		{
			if (p.x < pMin.x) pMin.x = p.x;
			if (p.y < pMin.y) pMin.y = p.y;
			if (p.x > pMax.x) pMax.x = p.x;
			if (p.y > pMax.y) pMax.y = p.y;
		}

		void Clear()
		{
			pMin.x = pMin.y = INF; 
			pMax.x = pMax.y = -INF;
		}
	};

	struct Grid2D
	{
		Grid2D(double _dx, double _dy);
		Grid2D(Grid2D &grid);
		~Grid2D();

		int dimx, dimy;
		double dx, dy;

		CellType GetType(int x, int y);
		CondData2D GetData(int x, int y);
		void Prepare(int frame, double substep);

		int LoadFromFile(char *filename);
		void TestPrint();

		BBox2D bbox;
	protected:
		//---------------------------- Borders ---------------------------------

		Point2D** points;
		Vec2D* velocities;
		int num_points;
		int num_frames;
		Vec2D* ComputeBorderVelocities(int frame, double substep);
		Point2D* ComputeSubBorder(int frame, double substep);
		//----------------------------------------------------------------------

		int *typeData;
		CondData2D *initData;

		void Init();

		void RasterLine(Point2D p1, Point2D p2, Vec2D v1, Vec2D v2, CellType color);
		void FloodFill(CellType color);

		void FillTestInitData(Vec2D startVel);
		
		void BuildBBox(int num_points, int num_frames, Point2D** points);
		void Build(int num_points, Point2D *points, Vec2D* vels);
		
		void SetType(int x, int y, CellType t);
		void SetData(int x, int y, CondData2D d);
	};
}