#pragma once

#include <stdio.h>

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
	};

	struct Grid3D
	{
		int dimx, dimy, dimz;
		double dx, dy, dz;

		Grid3D(double _dx, double _dy, double _dz);

		bool LoadFromFile(char *filename);

	protected:
		Node *nodes;
	};
}