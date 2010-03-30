#pragma once

#include "..\Common\Common.h"
#include "..\FluidSolver2D\Grid2D.h"

using namespace Common;

namespace FluidSolver3D
{
	enum NodeType { IN, OUT, BOUND };
	enum BCtype { NOSLIP, FREE };

	struct Node
	{
		NodeType type;

		BCtype bc;
		Vec3D v;
		double T;

		void SetBound(BCtype _bc, Vec3D _v, double _T)
		{
			type = BOUND;
			bc = _bc;
			v = _v;
			T = _T;
		}
	};

	struct Grid3D
	{
		int dimx, dimy, dimz;
		double dx, dy, dz;

		double baseT;

		Grid3D(double _dx, double _dy, double _dz, double _depth, double _baseT);
		~Grid3D();

		NodeType GetType(int i, int j, int k);
		Vec3D GetVel(int i, int j, int k);
		double GetT(int i, int j, int k);
		double GetFrameTime();
		
		void SetNodeVel(int i, int j, int k, Vec3D new_v);
		
		bool LoadFromFile(char *filename);
		void Prepare(double time);

	protected:
		Node*	nodes;		// all grid nodes

		FluidSolver2D::Grid2D *grid2D;		// 2D helper grid for borders
		double depth;						// depth
	};
}