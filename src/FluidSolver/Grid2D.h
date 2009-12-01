#pragma once

#include "Common.h"

#include <stdlib.h>
#include <math.h>

#include <vector>
#include <string>

#define INF				1e10
#define MAX_STR_SIZE	255
#define NUM_STEPS		40
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
		~Grid2D();

		int LoadFromFile(char *filename);
		void TestPrint();
	
	private:
		double dx, dy;
		int dimx, dimy;
		int *data;

		BBox2D bbox;

		void RasterLine(Point2D p1, Point2D p2, int steps, CellType color);
		void FloodFill(CellType color);

		void BuildBBox(int num_points, Point2D *points);
		void Build(int num_points, Point2D *points);

		CellType GetCell(int x, int y);
		void SetCell(int x, int y, CellType t);
	};
}