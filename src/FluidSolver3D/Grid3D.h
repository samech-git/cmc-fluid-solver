#pragma once

#include "..\Common\Structures.h"
#include "..\Common\IO.h"
#include "..\Common\Algorithms.h"

#include "..\FluidSolver2D\Grid2D.h"

using namespace Common;

namespace FluidSolver3D
{
	enum NodeType { 
		NODE_IN, 
		NODE_OUT, 
		NODE_BOUND, 
		NODE_VALVE 
	};
	
	enum BCtype { 
		BC_NOSLIP, 
		BC_FREE 
	};

	struct Node
	{
		NodeType type;

		BCtype bc_vel, bc_temp;
		Vec3D v;
		double T;

		void SetBound(BCtype _bc_vel, BCtype _bc_temp, Vec3D _v, double _T)
		{
			type = NODE_BOUND;
			bc_vel = _bc_vel;
			bc_temp = _bc_temp;
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
		BCtype GetBC_vel(int i, int j, int k);
		BCtype GetBC_temp(int i, int j, int k);
		Vec3D GetVel(int i, int j, int k);
		double GetT(int i, int j, int k);

		void SetNodeVel(int i, int j, int k, Vec3D new_v);

		double GetFrameTime();
		FluidSolver2D::Grid2D *GetGrid2D();
		
		bool LoadFromFile(char *filename);
		void Prepare(double time);

		void Grid3D::TestPrint(char *filename);

	protected:
		Node*	nodes;		// all grid nodes

		FluidSolver2D::Grid2D *grid2D;		// 2D helper grid for borders
		double depth;						// depth
	};
}