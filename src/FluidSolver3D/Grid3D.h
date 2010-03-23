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
	};

	struct BoundNode : Node
	{
		BCtype bc;
		Vec3D v;
		double T;
	};

	struct Grid3D
	{
		int dimx, dimy, dimz;
		double dx, dy, dz;

		Grid3D(double _dx, double _dy, double _dz, double _depth);
		~Grid3D();

		NodeType GetType(int i, int j, int k);

		bool LoadFromFile(char *filename);
		void Prepare(double time);

	protected:
		Node *nodes;		// all grid nodes

		FluidSolver2D::Grid2D *grid2D;		// 2D helper grid for borders
		double depth;						// depth
	};
}