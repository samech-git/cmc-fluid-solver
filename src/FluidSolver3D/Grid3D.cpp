#include "Grid3D.h"

namespace FluidSolver3D
{
	Grid3D::Grid3D(double _dx, double _dy, double _dz, double _depth, double _baseT) : 
		dx(_dx), dy(_dy), dz(_dz), depth(_depth), baseT(_baseT), nodes(NULL)
	{
		grid2D = new FluidSolver2D::Grid2D(dx, dy, 0, true);
	}

	Grid3D::~Grid3D()
	{
		if (nodes != NULL) delete [] nodes;
		if (grid2D != NULL) delete grid2D;
	}

	NodeType Grid3D::GetType(int i, int j, int k)
	{
		return nodes[i * dimy * dimz + j * dimz + k].type;
	}

	Vec3D Grid3D::GetVel(int i, int j, int k)
	{
		return nodes[i * dimy * dimz + j * dimz + k].v;
	}

	double Grid3D::GetT(int i, int j, int k)
	{
		return nodes[i * dimy * dimz + j * dimz + k].T;
	}

	double Grid3D::GetFrameTime()
	{
		int frames = grid2D->GetFramesNum();
		double length = grid2D->GetCycleLenght();
		return (length / frames);
	}

	void Grid3D::SetNodeVel(int i, int j, int k, Vec3D new_v)
	{
		nodes[i * dimy * dimz + j * dimz + k].v = new_v;
	}

	bool Grid3D::LoadFromFile(char *filename)
	{
		if (grid2D->LoadFromFile(filename))
		{
			dimx = grid2D->dimx;
			dimy = grid2D->dimy;
			dimz = (int)ceil(depth / dz) + 1;
			nodes = new Node[dimx * dimy * dimz];
			return true;
		}
		else
			return false;
	}

	void Grid3D::Prepare(double time)
	{
		grid2D->Prepare(time);

		memset(nodes, 0, sizeof(Node) * dimx * dimy * dimz);
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
			{
				int ind;
				if (grid2D->GetType(i, j) == FluidSolver2D::OUT)
				{
					for (int k = 0; k < dimz; k++)
					{
						ind = i * dimy * dimz + j * dimz + k;
						nodes[ind].type = OUT;
					}
				}
				else
				{
					// set up & bottom bounds
					nodes[i * dimy * dimz + j * dimz + 0].type = OUT;
					nodes[i * dimy * dimz + j * dimz + dimz-1].type = OUT;

					nodes[i * dimy * dimz + j * dimz + 1].SetBound(NOSLIP, Vec3D(0.0, 0.0, 0.0), 0.0);
					nodes[i * dimy * dimz + j * dimz + dimz-2].SetBound(NOSLIP, Vec3D(0.0, 0.0, 0.0), 0.0);
					
					for (int k = 2; k < dimz-2; k++)
					{
						ind = i * dimy * dimz + j * dimz + k;
						switch (grid2D->GetType(i, j))
						{
						case FluidSolver2D::BOUND:
						case FluidSolver2D::VALVE:
							nodes[ind].SetBound(NOSLIP, Vec3D(grid2D->GetData(i, j).vel.x, grid2D->GetData(i, j).vel.y, 0.0), grid2D->GetData(i, j).T);
							break;
						case FluidSolver2D::IN:
							nodes[ind].type = IN;
							break;
						}
					}
				}
			}
	}
}