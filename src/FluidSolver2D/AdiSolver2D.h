#pragma once

#include "Solver2D.h"

#define ERR_THRESHOLD		0.1
#define MAX_GLOBAL_ITERS	100

namespace FluidSolver2D
{
	enum VarType { type_U, type_V, type_T };
	enum DirType { X, Y };

	struct Segment2D
	{
		int posx, posy;
		int endx, endy;
		int size;
		DirType dir; 
	};

	class AdiSolver2D : public Solver2D
	{
	public:
		void Init(Grid2D* _grid, FluidParams &_params);
		void TimeStep(double dt, int num_global, int num_local);

		AdiSolver2D();
		~AdiSolver2D();
	
	private:
		vector<Segment2D> listX, listY;

		TimeLayer2D *temp, *half, *next_local;
		double *a, *b, *c, *d, *x;

		void SolveTridiagonal(double *a, double *b, double *c, double *d, double *x, int num);

		void BuildMatrix(double dt, int i, int j, VarType var, DirType dir, double *a, double *b, double *c, double *d, int n, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *temp_local);
		void ApplyBC0(int i, int j, VarType var, double &b0, double &c0, double &d0);
		void ApplyBC1(int i, int j, VarType var, double &a1, double &b1, double &d1);
		
		void CreateSegments();	
		void SolveSegment(double dt, Segment2D seg, VarType var, DirType dir, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *temp_local, TimeLayer2D *next_local);
		void UpdateSegment(double *x, Segment2D seg, VarType var, TimeLayer2D *layer);
		
		void SolveDirection(double dt, int num_local, vector<Segment2D> &list, TimeLayer2D *cur, TimeLayer2D *temp, TimeLayer2D *next);

		void FreeMemory();
	};
}