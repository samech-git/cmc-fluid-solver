#include "Grid2D.h"

namespace FluidSolver
{
	Grid2D::Grid2D(double _dx, double _dy) : dx(_dx), dy(_dy), curData(NULL), nextData(NULL) {	}

	Grid2D::Grid2D(Grid2D &grid) : dx(grid.dx), dy(grid.dy), curData(NULL)
	{
		dimx = grid.dimx;
		dimy = grid.dimy;

		curData = new CondData2D[dimx * dimy];
		memcpy(curData, grid.curData, dimx * dimy * sizeof(CondData2D));

		nextData = new CondData2D[dimx * dimy];
		memcpy(nextData, grid.nextData, dimx * dimy * sizeof(Vec2D));
	}

	Grid2D::~Grid2D()
	{
		if (curData != NULL) delete [] curData;
		if (nextData != NULL) delete [] nextData;

		for (int i=0; i<num_frames; i++)
			frames[i].Dispose();
		delete[] frames;
	}

	inline CellType Grid2D::GetType(int x, int y)
	{
		return curData[x * dimy + y].cell;
	}

	CondData2D Grid2D::GetData(int x, int y)
	{
		return curData[x * dimy + y];
	}

	inline void Grid2D::SetType(int x, int y, CellType t)
	{
		curData[x * dimy + y].cell = t;
	}

	inline void Grid2D::SetData(int x, int y, CondData2D d)
	{
		curData[x * dimy + y] = d;
	}

	void Grid2D::SetFieldData(int x, int y, CondData2D d)
	{
		nextData[x * dimy + y] = d;
	}

	void ReadPoint2D(FILE *file, Point2D &p)
	{
		string str = "";
		char c;	

		// read line
		fscanf_s(file, "%c", &c);
		while (c == '\n' || c == ' ') fscanf_s(file, "%c", &c);
		while (c != '\n')
		{
			str += c;	
			fscanf_s(file, "%c", &c);
		}

		// replace ',' with '.' if necessary
		string::size_type pos = 0;
		string::size_type found;
		while ((found = str.find(',', pos)) != string::npos)
		{
			str.replace(found, 1, 1, '.');
			pos = found;
		}

		// separate 2 values
		pos = 0;
		found = str.find(' ', pos);
		string s1 = str.substr(0, found);
		string s2 = str.substr(found+1, string::npos);

		// convert to doubles
		p.x = atof(s1.c_str());
		p.y = atof(s2.c_str());
	}

	void Grid2D::BuildBBox(int num_frames, FrameInfo* frames)
	{
		bbox.Clear();
		for (int j = 0; j < num_frames; j++)
			for (int i = 0; i < frames[j].NumShapes; i++)
				for (int k = 0; k < frames[j].Shapes[i].NumPoints; k++)
					bbox.AddPoint(frames[j].Shapes[i].Points[k]);

		double wx = bbox.pMax.x - bbox.pMin.x;
		double wy = bbox.pMax.y - bbox.pMin.y;

		bbox.pMin.x -= wx * 0.02;
		bbox.pMin.y -= wy * 0.02;
		bbox.pMax.x += wx * 0.02;
		bbox.pMax.y += wy * 0.02;

#ifdef _DEBUG
		printf("bbox x range: [%f, %f]\n", bbox.pMin.x, bbox.pMax.x);
		printf("bbox y range: [%f, %f]\n", bbox.pMin.y, bbox.pMax.y);
#endif
	}

	VecTN Grid2D::GetTangentNormal(Vec2D vector, Vec2D orientation)
	{
		double l = (vector.x * orientation.x + vector.y * orientation.y) / (orientation.x * orientation.x + orientation.y * orientation.y);
		Vec2D t(orientation.x * l, orientation.y * l);
		Vec2D n(vector.x - t.x, vector.y - t.y);
		return VecTN(t, n);
	}

#define PROCESS(ij) {if (nextData[ij].cell != OUT) { v.x += nextData[ij].vel.x; v.y += nextData[ij].vel.y; k++; }}

	Vec2D Grid2D::GetBounfVelocitie(int x, int y)
	{
		int ij = x * dimy + y;
		Vec2D v(0, 0);
		int k = 0;

		PROCESS(ij - dimy-1);
		PROCESS(ij - dimy);
		PROCESS(ij - dimy+1);
		PROCESS(ij - 1);
		PROCESS(ij);
		PROCESS(ij + 1);
		PROCESS(ij + dimy-1);
		PROCESS(ij + dimy);
		PROCESS(ij + dimy+1);

		if (k != 0)
		{
			v.x /= k;
			v.y /= k;
		}
		return v;
	}

	void Grid2D::RasterLine(Point2D p1, Point2D p2, Vec2D v1, Vec2D v2, CellType color)
    {
		Vec2D orientation(p2.x - p1.x, p2.y - p1.y);
		int steps = (int)max(abs(orientation.x), abs(orientation.y)) + 1;
        Point2D dp((orientation.x) / steps, (orientation.y) / steps);		// divide a segment into parts
		Vec2D dv((v2.x - v1.x) / steps, (v2.y - v1.y) / steps);		// divide a segment into parts
		
		Point2D p = p1;
		Vec2D v = v1;

		// go through the whole segment
        for (int i = 0; i <= steps; i++) 
        {
            int x = (int)p.x;
            int y = (int)p.y;
			
			Vec2D bv = GetBounfVelocitie(x, y);
			VecTN vtn = GetTangentNormal(v, orientation);
			VecTN btn = GetTangentNormal(bv, orientation);

#ifndef BC_NOSLIP
			SetData(x, y, CondData2D(NOSLIP, color, Vec2D(vtn.normal.x + btn.tangent.x, vtn.normal.y + btn.tangent.y), 1.0));
#else
			SetData(x, y, CondData2D(NOSLIP, color, Vec2D(v.x, v.y), 1.0));
#endif

			p.x += dp.x;
            p.y += dp.y;
			v.x += dv.x;
            v.y += dv.y;
        }
    }

	void Grid2D::FloodFill(CellType color)
    {
		// go through all the rows
        for (int j = 0; j < dimy; j++) 
        {
			// left to right
            int i = 0;
            while (i < dimx && GetType(i, j) == IN) 
            {
                SetType(i, j, color);
                i++;
            }

			// right to left
            i = dimx - 1;
            while (GetType(i, j) == IN && i >= 0) 
            {
                SetType(i, j, color);
                i--;
            }
        }

		// go through all the column
		for (int i = 0; i < dimx; i++) 
		{
			// up to down
			int j = 0;
			while (j < dimy && GetType(i, j) == IN) 
			{
				SetType(i, j, color);
				j++;
			}

			// down to up
			j = dimy - 1;
			while (GetType(i, j) == IN && j >= 0) 
			{
				SetType(i, j, color);
				j--;
			}
		}
    }

	void Grid2D::Init()
	{
		BuildBBox(num_frames, frames);
	
		dimx = (int)ceil((bbox.pMax.x - bbox.pMin.x) / dx) + 1;
		dimy = (int)ceil((bbox.pMax.y - bbox.pMin.y) / dy) + 1;

		// allocate data
		int size = dimx * dimy;
		curData = new CondData2D[size];
		nextData = new CondData2D[size];

		for (int i=0; i<size; i++)
		{
			nextData[i].T = 0;
			nextData[i].type = NONE;
			nextData[i].cell = OUT;
			nextData[i].vel.x = 0;
			nextData[i].vel.y = 0;
		}



		// convert physical coordinates to the grid coordinates
		for (int j = 0; j < num_frames; j++)
			for (int i = 0; i < frames[j].NumShapes; i++)
				for (int k = 0; k < frames[j].Shapes[i].NumPoints; k++)
				{
					frames[j].Shapes[i].Points[k].x = (frames[j].Shapes[i].Points[k].x - bbox.pMin.x) / dx;
					frames[j].Shapes[i].Points[k].y = (frames[j].Shapes[i].Points[k].y - bbox.pMin.y) / dy;
				}

	}

	void Grid2D::Build(FrameInfo frame)
	{
		// mark all cells as inner 
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				SetType(i, j, IN);
     
		// rasterize lines
		for (int j=0; j<frame.NumShapes; j++)
			for (int i = 0; i < frame.Shapes[j].NumPoints - 1; i++) 
				RasterLine(frame.Shapes[j].Points[i], frame.Shapes[j].Points[i+1],
						   frame.Shapes[j].Velocities[i], frame.Shapes[j].Velocities[i+1],
						   BOUND);
        FloodFill(OUT); 

		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
			{
				CellType c = GetType(i, j);
				switch (c)
				{
					case IN: case OUT: SetData(i, j, CondData2D(NONE, c, Vec2D(0, 0), 1.0)); break; 
				}
			}
	}

	int Grid2D::LoadFromFile(char *filename)
	{
		FILE *file = NULL;
		if (fopen_s(&file, filename, "r") != 0)
		{
			printf("Error: cannot open file \"%s\" \n", filename);
			return ERR;
		}

		fscanf_s(file, "%i", &num_frames);	// currently not used
		Point2D p;
		int temp;
		frames = new FrameInfo[num_frames];
		
		for (int j=0; j<num_frames; j++)
		{
			
			fscanf_s(file, "%f", &(frames[j].Duration));
			frames[j].Duration = 0.035;
			fscanf_s(file, "%i", &temp);
			frames[j].Init(temp);
			for (int i = 0; i<frames[j].NumShapes; i++)
			{
				fscanf_s(file, "%i", &temp);
				frames[j].Shapes[i].Init(temp);
				for (int k=0; k<frames[j].Shapes[i].NumPoints; k++)
				{
					ReadPoint2D(file, p);
					frames[j].Shapes[i].Points[k].x = p.x / 1000;
					frames[j].Shapes[i].Points[k].y = p.y / 1000;
				}

				char str[8];
				fscanf_s(file, "%s", str, 8);

				p.x = p.y = 0;
				if (str[0] == 'M')
				{
					frames[j].Shapes[i].Active = true;
					ReadPoint2D(file, p);
				}
				else
					frames[j].Shapes[i].Active = false;
				for (int k=0; k<frames[j].Shapes[i].NumPoints; k++)
				{
					frames[j].Shapes[i].Velocities[k].x = p.x / 1000;
					frames[j].Shapes[i].Velocities[k].y = p.y / 1000;
				}
			}
		}
	
		for (int j=0; j<num_frames; j++)
			ComputeBorderVelocities(j);
		
		Init();
		return OK;
	}


	void Grid2D::ComputeBorderVelocities(int frame)
	{
		int nextframe = (frame + 1) % num_frames;
		double m = 1 / frames[frame].Duration;
		for (int i=0; i<frames[frame].NumShapes; i++)
			if (!frames[frame].Shapes[i].Active)
				for (int k=0; k<frames[frame].Shapes[i].NumPoints; k++)
				{
					frames[nextframe].Shapes[i].Velocities[k].x = 
						(frames[nextframe].Shapes[i].Points[k].x - frames[frame].Shapes[i].Points[k].x) * m;
					frames[nextframe].Shapes[i].Velocities[k].y = 
						(frames[nextframe].Shapes[i].Points[k].y - frames[frame].Shapes[i].Points[k].y) * m;
				}
			else
				for (int k=0; k<frames[frame].Shapes[i].NumPoints; k++)
				{
					frames[nextframe].Shapes[i].Velocities[k].x += 
						(frames[frame].Shapes[i].Points[k].x - frames[nextframe].Shapes[i].Points[k].x) * m;
					frames[nextframe].Shapes[i].Velocities[k].y += 
						(frames[frame].Shapes[i].Points[k].y - frames[nextframe].Shapes[i].Points[k].y) * m;
				}
	}


	FrameInfo Grid2D::ComputeSubframe(int frame, double substep)
	{
		int framep1 = (frame + 1) % num_frames;

		FrameInfo res;
		res.Duration = 0;

		res.Init(frames[frame].NumShapes);
		for (int i=0; i<res.NumShapes; i++)
		{
			res.Shapes[i].Init(frames[frame].Shapes[i].NumPoints);
			for (int k=0; k<res.Shapes[i].NumPoints; k++)
			{
				res.Shapes[i].Points[k].x = frames[frame].Shapes[i].Points[k].x * (1 - substep) + 
											frames[framep1].Shapes[i].Points[k].x * substep;
				res.Shapes[i].Points[k].y = frames[frame].Shapes[i].Points[k].y * (1 - substep) + 
											frames[framep1].Shapes[i].Points[k].y * substep;

				res.Shapes[i].Velocities[k].x = frames[frame].Shapes[i].Velocities[k].x * (1 - substep) + 
												frames[framep1].Shapes[i].Velocities[k].x * substep;
				res.Shapes[i].Velocities[k].y = frames[frame].Shapes[i].Velocities[k].y * (1 - substep) + 
												frames[framep1].Shapes[i].Velocities[k].y * substep;
			}
		}
		return res;
	}

	void Grid2D::Prepare(int frame, double substep)
	{
		FrameInfo tempframe = ComputeSubframe(frame % num_frames, substep);
		Build(tempframe);
		tempframe.Dispose();
	}

	void Grid2D::Prepare(double time)
	{
		double* a = new double[num_frames + 1];
		a[0] = 0;
		for (int i=1; i<=num_frames; i++)
			a[i] = a[i-1] + frames[i-1].Duration;

		double r_time = fmod(time, a[num_frames]);
		int frame = 0;
		for (int i=1; i<num_frames; i++)
			if (a[i] < r_time) frame = i;
		double substep = (r_time - a[frame]) / (a[frame + 1] - a[frame]);
		delete[] a;

		Prepare(frame, substep);
	}


	void Grid2D::TestPrint()
	{
		printf("grid view:\n");
		printf("%i %i\n", dimx, dimy);
		for (int i = 0; i < dimx; i++)
		{
			for (int j = 0; j < dimy; j++)
			{
				CellType t = GetType(i, j);
				switch (t)
				{
				case IN: printf(" "); break;
				case OUT: printf("."); break;
				case BOUND: printf("#"); break;
				case VALVE: printf("+"); break;
				}
			}
			printf("\n");
		}
	}
}