#pragma once

#include "Solver3D.h"

#define ERR_THRESHOLD		0.1
#define	MAX_GLOBAL_ITERS	10

namespace FluidSolver3D
{
	enum VarType { type_U, type_V, type_W, type_T };
	enum DirType { X, Y, Z };

	struct Segment3D
	{
		int posx, posy, posz;
		int endx, endy, endz;
		int size;
		DirType dir; 
	};

	class AdiSolver3D : public Solver3D
	{
	public:
		AdiSolver3D();
		~AdiSolver3D();

		void Init(Grid3D* _grid, FluidParams &_params);
		void TimeStep(double dt, int num_global, int num_local);

	private:
		vector<Segment3D> listX, listY, listZ;

		TimeLayer3D *temp, *half1, *half2, *next_local;
		double *a, *b, *c, *d, *x;

		void BuildMatrix(double dt, int i, int j, int k, VarType var, DirType dir, double *a, double *b, double *c, double *d, int n, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *temp_local);
		void ApplyBC0(int i, int j, int k, VarType var, double &b0, double &c0, double &d0);
		void ApplyBC1(int i, int j, int k, VarType var, double &a1, double &b1, double &d1);
		
		void CreateSegments();
		void SolveSegment(double dt, Segment3D seg, VarType var, DirType dir, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *temp_local, TimeLayer3D *next_local);
		void UpdateSegment(double *x, Segment3D seg, VarType var, TimeLayer3D *layer);
		
		void SolveDirection(double dt, int num_local, vector<Segment3D> &list, TimeLayer3D *cur, TimeLayer3D *temp, TimeLayer3D *next);

		void FreeMemory();
	};
}