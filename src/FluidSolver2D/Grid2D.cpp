#include "Grid2D.h"
#include "IO.h"

namespace FluidSolver2D
{
	Grid2D::Grid2D(double _dx, double _dy, double _startT, bool _bc_noslip, double _bc_strength) : dx(_dx), dy(_dy), startT(_startT), bc_noslip(_bc_noslip), bc_strength(_bc_strength), curData(NULL), nextData(NULL) {	}

	Grid2D::Grid2D(Grid2D &grid) : dx(grid.dx), dy(grid.dy), startT(grid.startT), bc_noslip(grid.bc_noslip), bc_strength(grid.bc_strength), curData(NULL)
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

	int ExtractInt(char* str)
	{
		int result = 0;
		for(int i=0; i>=0; i++)
		{
			char c = str[i];
			if (c >= '0' && c<= '9')
				result = result*10 + c - '0';
			if (c == 0)
				return result;
		}
		return result;
	}

	VecTN Grid2D::GetTangentNormal(Vec2D vector, Vec2D orientation)
	{
		double l = (vector.x * orientation.x + vector.y * orientation.y) / (orientation.x * orientation.x + orientation.y * orientation.y);
		Vec2D t(orientation.x * l, orientation.y * l);
		Vec2D n(vector.x - t.x, vector.y - t.y);
		return VecTN(t, n);
	}

#define PROCESS(ij) {if (nextData[ij].cell != CELL_OUT) { v.x += nextData[ij].vel.x; v.y += nextData[ij].vel.y; k++; }}

	Vec2D Grid2D::GetBounfVelocity(int x, int y)
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
			
			Vec2D bv = GetBounfVelocity(x, y);
			VecTN vtn = GetTangentNormal(v, orientation);
			VecTN btn = GetTangentNormal(bv, orientation);

			double l = sqrt(v.x*v.x + v.y*v.y);

			//bc_strength;
			if (bc_noslip)
				SetData(x, y, CondData2D(COND_NOSLIP, color, Vec2D(v.x, v.y), startT));
			else
				SetData(x, y, CondData2D(COND_NOSLIP, color, 
						Vec2D(vtn.normal.x + (btn.tangent.x*bc_strength + vtn.tangent.x*(1-bc_strength)), 
							  vtn.normal.y + (btn.tangent.y*bc_strength + vtn.tangent.y*(1-bc_strength))), 
						startT));

			p.x += dp.x;
            p.y += dp.y;
			v.x += dv.x;
            v.y += dv.y;
        }
    }

	void Grid2D::RasterField(Field2D field)
	{
		for (int j=1; j<dimy-1; j++)
			for (int i=1; i<dimx-1; i++)
			{
				double x = bbox.pMin.x + i*dx;
				double y = bbox.pMin.y + j*dy;
				Vec2D v = field.GetVelocity(x, y);
				if (v.x != 0 || v.y != 0)
					SetData(i, j, CondData2D(COND_NOSLIP, CELL_BOUND, v, startT));
			};
	}

	void Grid2D::FloodFill(CellType color)
    {
		int *queue = new int[dimx * dimy * 2];
		int cur = -1;

		const int neighborPos[4][2] = { {-1, 0}, {1, 0}, {0, -1}, {0, 1} };

		// we know that this cell is of our color type
		int last = 0;
		queue[0] = 0;
		queue[1] = 0;
		SetType(0, 0, color);

		// run wave
		while (cur < last)
		{
			// get next index
			cur++;
			int i = queue[cur * 2 + 0];
			int j = queue[cur * 2 + 1];

			// add neighbours
			for (int k = 0; k < 4; k++)
			{
				int next_i = i + neighborPos[k][0];
				int next_j = j + neighborPos[k][1];

				if ((next_i >= 0) && (next_i < dimx) && (next_j >= 0) && (next_j < dimy))
					if (GetType(next_i, next_j) == CELL_IN) 
					{
						last++;
						queue[last * 2 + 0] = next_i;
						queue[last * 2 + 1] = next_j;
						SetType(next_i, next_j, color);
					}
			}
		}

		delete [] queue;
	}


	void Grid2D::Init()
	{
		bbox.Build(num_frames, frames);
	
		dimx = (int)ceil((bbox.pMax.x - bbox.pMin.x) / dx) + 1;
		dimy = (int)ceil((bbox.pMax.y - bbox.pMin.y) / dy) + 1;

		// allocate data
		int size = dimx * dimy;
		curData = new CondData2D[size];
		nextData = new CondData2D[size];

		for (int i=0; i<size; i++)
		{
			nextData[i].T = 0;
			nextData[i].type = COND_NONE;
			nextData[i].cell = CELL_OUT;
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

	void Grid2D::Build(FrameInfo2D frame)
	{
		// mark all cells as inner 
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				SetType(i, j, CELL_IN);
     
		// rasterize lines (VALVE)
		for (int j=0; j<frame.NumShapes; j++)
			for (int i = 0; i < frame.Shapes[j].NumPoints - 1; i++) 
				if( frame.Shapes[j].Active )
					RasterLine(frame.Shapes[j].Points[i], frame.Shapes[j].Points[i+1],
							   frame.Shapes[j].Velocities[i], frame.Shapes[j].Velocities[i+1],
							   CELL_VALVE);
		
		// rasterize lines (BOUND)
		for (int j=0; j<frame.NumShapes; j++)
			for (int i = 0; i < frame.Shapes[j].NumPoints - 1; i++) 
				if( !frame.Shapes[j].Active )
					RasterLine(frame.Shapes[j].Points[i], frame.Shapes[j].Points[i+1],
							   frame.Shapes[j].Velocities[i], frame.Shapes[j].Velocities[i+1],
							   CELL_BOUND);
			
        FloodFill(CELL_OUT); 
		RasterField(frame.Field);

		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
			{
				CellType c = GetType(i, j);
				switch (c)
				{
					case CELL_IN: case CELL_OUT: SetData(i, j, CondData2D(COND_NONE, c, Vec2D(0, 0), startT)); break; 
				}
			}
	}

	bool Grid2D::LoadFromFile(char *filename, char *fieldname)
	{
		FILE *file = NULL;
		if (fopen_s(&file, filename, "r") != 0)
		{
			printf("Error: cannot open file \"%s\" \n", filename);
			return false;
		}

		fscanf_s(file, "%i", &num_frames);	// currently not used
		Point2D p;
		int temp;
		frames = new FrameInfo2D[num_frames];
		
		for (int j=0; j<num_frames; j++)
		{
			float dur;
			fscanf_s(file, "%f", &dur);
			frames[j].Duration = dur;
			fscanf_s(file, "%i", &temp);
			frames[j].Init(temp);
			for (int i = 0; i<frames[j].NumShapes; i++)
			{
				fscanf_s(file, "%i", &temp);
				frames[j].Shapes[i].Init(temp);
				for (int k=0; k<frames[j].Shapes[i].NumPoints; k++)
				{
					ReadPoint2D(file, p);
					frames[j].Shapes[i].Points[k].x = p.x * GRID_SCALE_FACTOR;
					frames[j].Shapes[i].Points[k].y = p.y * GRID_SCALE_FACTOR;
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
					frames[j].Shapes[i].Velocities[k].x = p.x * GRID_SCALE_FACTOR;
					frames[j].Shapes[i].Velocities[k].y = p.y * GRID_SCALE_FACTOR;
				}
			}
		}
		fclose(file);

		if (strcmp(fieldname, "") != 0)
		{
			if (fopen_s(&file, fieldname, "r") != 0)
			{
				printf("Error: cannot open file \"%s\" \n", filename);
				return false;
			}

			float minx, miny;
			float maxx, maxy;
			float dx, dy;
			int nx, ny;
			
			fscanf_s(file, "%f %f %f %f\n", &minx, &miny, &maxx, &maxy);
			fscanf_s(file, "%f %f %i %i\n", &dx, &dy, &nx, &ny);

			minx *= GRID_SCALE_FACTOR;
			miny *= GRID_SCALE_FACTOR;
			maxx *= GRID_SCALE_FACTOR;
			maxy *= GRID_SCALE_FACTOR;

			dx *= GRID_SCALE_FACTOR;
			dy *= GRID_SCALE_FACTOR;
				
			for(;;)
			{
				char str[255];
				ReadLine(file, str, 255);
				if (str[0] != 'F') break;

				int frame = ExtractInt(str);
				frames[frame].Field.Init(minx, miny, dx, dy, nx, ny);
				ReadLine(file, str, 255);

				float x, y;
				for (int j=0; j<ny; j++)
					for (int i=0; i<nx; i++)
					{
						fscanf_s(file, "%f %f", &x, &y);
						frames[frame].Field.Data[j*nx + i] = Vec2D(x, y);
					}
				ReadLine(file, str, 255);
			}
			fclose(file);
		}
	
		for (int j=0; j<num_frames; j++)
			ComputeBorderVelocities(j);
		
		Init();
		return true;
	}


	void Grid2D::ComputeBorderVelocities(int frame)
	{
		int nextframe = (frame + 1) % num_frames;
		double m = 1 / frames[frame].Duration;
		
		//double maxl = 0;//!!!!!!

		for (int i=0; i<frames[frame].NumShapes; i++)
			if (!frames[frame].Shapes[i].Active)
				for (int k=0; k<frames[frame].Shapes[i].NumPoints; k++)
				{
					frames[nextframe].Shapes[i].Velocities[k].x = 
						(frames[nextframe].Shapes[i].Points[k].x - frames[frame].Shapes[i].Points[k].x) * m;
					frames[nextframe].Shapes[i].Velocities[k].y = 
						(frames[nextframe].Shapes[i].Points[k].y - frames[frame].Shapes[i].Points[k].y) * m;

					//double x = frames[nextframe].Shapes[i].Velocities[k].x;
					//double y = frames[nextframe].Shapes[i].Velocities[k].y;
					//double l = sqrt(x*x + y*y);
					//if (l > maxl) maxl = l;//!!!!!!
				}
			else
				for (int k=0; k<frames[frame].Shapes[i].NumPoints; k++)
				{
					frames[nextframe].Shapes[i].Velocities[k].x += 
						(frames[frame].Shapes[i].Points[k].x - frames[nextframe].Shapes[i].Points[k].x) * m;
					frames[nextframe].Shapes[i].Velocities[k].y += 
						(frames[frame].Shapes[i].Points[k].y - frames[nextframe].Shapes[i].Points[k].y) * m;

					//double x = frames[nextframe].Shapes[i].Velocities[k].x;
					//double y = frames[nextframe].Shapes[i].Velocities[k].y;
					//double l = sqrt(x*x + y*y);
					//if (l > maxl) maxl = l;//!!!!!!
				}
	}


	FrameInfo2D Grid2D::ComputeSubframe(int frame, double substep)
	{
		int framep1 = (frame + 1) % num_frames;

		FrameInfo2D res;
		res.Duration = 0;

		double isubstep = 1 - substep;

		res.Init(frames[frame].NumShapes);
		for (int i=0; i<res.NumShapes; i++)
		{
			res.Shapes[i].Init(frames[frame].Shapes[i].NumPoints);
			for (int k=0; k<res.Shapes[i].NumPoints; k++)
			{
				res.Shapes[i].Points[k].x = frames[frame].Shapes[i].Points[k].x * isubstep + 
											frames[framep1].Shapes[i].Points[k].x * substep;
				res.Shapes[i].Points[k].y = frames[frame].Shapes[i].Points[k].y * isubstep + 
											frames[framep1].Shapes[i].Points[k].y * substep;

				res.Shapes[i].Velocities[k].x = frames[frame].Shapes[i].Velocities[k].x * isubstep + 
												frames[framep1].Shapes[i].Velocities[k].x * substep;
				res.Shapes[i].Velocities[k].y = frames[frame].Shapes[i].Velocities[k].y * isubstep + 
												frames[framep1].Shapes[i].Velocities[k].y * substep;
				
				res.Shapes[i].Active = frames[frame].Shapes[i].Active;
			}
		}

		if (frames[frame].Field.Correlate(frames[framep1].Field))
		{
			res.Field.Init(frames[frame].Field);
			int nx = res.Field.Nx;
			int ny = res.Field.Ny;

			for (int j=0; j<ny; j++)
				for (int i=0; i<nx; i++)
				{
					int t = j*nx + i;

					double x1 = frames[frame].Field.Data[t].x;
					double y1 = frames[frame].Field.Data[t].y;
					double x2 = frames[framep1].Field.Data[t].x;
					double y2 = frames[framep1].Field.Data[t].y;
					double x = 0;
					double y = 0;

					if ( (x1 != 0 || y1 != 0) && (x2 != 0 || y2 != 0) )
					{
						x = x1 * isubstep + x2 * substep;
						y = y1 * isubstep + y2 * substep;
					}

					res.Field.Data[t] = Vec2D(x, y);
				};
		}

		return res;
	}

	void Grid2D::Prepare(int frame, double substep)
	{
		FrameInfo2D tempframe = ComputeSubframe(frame % num_frames, substep);
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

	double Grid2D::GetCycleLenght()
	{
		double res = 0;
		for (int i=0; i<num_frames; i++)
			res += frames[i].Duration;
		return res;
	}

	int Grid2D::GetFramesNum()
	{
		return num_frames;
	}

	int Grid2D::GetFrame(double time)
	{
		double* a = new double[num_frames + 1];
		a[0] = 0;
		for (int i=1; i<=num_frames; i++)
			a[i] = a[i-1] + frames[i-1].Duration;

		double r_time = fmod(time, a[num_frames]);
		int frame = 0;
		for (int i=1; i<num_frames; i++)
			if (a[i] < r_time) frame = i;

		return frame;
	}

	float Grid2D::GetLayerTime(double t)
	{
		double* a = new double[num_frames + 1];
		a[0] = 0;
		for (int i=1; i<=num_frames; i++)
			a[i] = a[i-1] + frames[i-1].Duration;

		double r_time = fmod(t, a[num_frames]);
		int frame = 0;
		for (int i=1; i<num_frames; i++)
			if (a[i] < r_time) frame = i;

		return (float)(a[frame+1] - r_time);
	}


	void Grid2D::TestPrint(char *filename)
	{
		FILE *file = NULL;
		fopen_s(&file, filename, "w");
		fprintf(file, "grid view:\n");
		fprintf(file, "%i %i\n", dimx, dimy);
		for (int i = 0; i < dimx; i++)
		{
			for (int j = 0; j < dimy; j++)
			{
				CellType t = GetType(i, j);
				switch (t)
				{
				case CELL_IN: fprintf(file, " "); break;
				case CELL_OUT: fprintf(file, "."); break;
				case CELL_BOUND: fprintf(file, "#"); break;
				case CELL_VALVE: fprintf(file, "+"); break;
				}
			}
			fprintf(file, "\n");
		}
		fclose(file);
	}
}