#include "Grid2D.h"

namespace FluidSolver
{
	Grid2D::Grid2D(double _dx, double _dy) : dx(_dx), dy(_dy), typeData(NULL), initData(NULL) {	}

	Grid2D::Grid2D(Grid2D &grid) : dx(grid.dx), dy(grid.dy), typeData(NULL), initData(NULL)
	{
		dimx = grid.dimx;
		dimy = grid.dimy;

		typeData = new int[dimx * dimy];
		memcpy(typeData, grid.typeData, dimx * dimy * sizeof(int));

		initData = new CondData2D[dimx * dimy];
		memcpy(initData, grid.initData, dimx * dimy * sizeof(CondData2D));
	}

	Grid2D::~Grid2D()
	{
		if (typeData != NULL) delete [] typeData;
		if (initData != NULL) delete [] initData;
		if (velocities != NULL) delete[] velocities;
		if (points != NULL)
		{
			for (int i=0; i<num_frames; i++)
				if (points[i] != NULL) delete[] points[i];
			delete[] points;
		}
	}

	inline CellType Grid2D::GetType(int x, int y)
	{
		return (CellType)typeData[x * dimy + y];
	}

	CondData2D Grid2D::GetData(int x, int y)
	{
		return initData[x * dimy + y];
	}

	inline void Grid2D::SetType(int x, int y, CellType t)
	{
		typeData[x * dimy + y] = t;
	}

	inline void Grid2D::SetData(int x, int y, CondData2D d)
	{
		initData[x * dimy + y] = d;
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

	void Grid2D::BuildBBox(int num_points, int num_frames, Point2D** points)
	{
		bbox.Clear();
		for (int j = 0; j < num_frames; j++)
			for (int i = 0; i < num_points; i++)
				bbox.AddPoint(points[j][i]);

		bbox.pMin.x -= BBOX_OFFSET;
		bbox.pMin.y -= BBOX_OFFSET;
		bbox.pMax.x += BBOX_OFFSET;
		bbox.pMax.y += BBOX_OFFSET;

#ifdef _DEBUG
		printf("bbox x range: [%f, %f]\n", bbox.pMin.x, bbox.pMax.x);
		printf("bbox y range: [%f, %f]\n", bbox.pMin.y, bbox.pMax.y);
#endif
	}

	void Grid2D::RasterLine(Point2D p1, Point2D p2, Vec2D v1, Vec2D v2, int steps, CellType color)
    {
        Point2D dp((p2.x - p1.x) / steps, (p2.y - p1.y) / steps);		// divide a segment into parts
		Vec2D dv((v2.x - v1.x) / steps, (v2.y - v1.y) / steps);		// divide a segment into parts
		
		Point2D startp = p1;
		Vec2D startv = v1;

		// go through the whole segment
        for (int i = 0; i <= steps; i++) 
        {
            int x = (int)startp.x;
            int y = (int)startp.y;

			SetType(x, y, color);
			SetData(x, y, CondData2D(NOSLIP, startv, 1.0));

			startp.x += dp.x;
            startp.y += dp.y;
			startv.x += dv.x;
            startv.y += dv.y;
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
    }

	void Grid2D::Init()
	{
		BuildBBox(num_points, num_frames, points);
	
		dimx = (int)ceil((bbox.pMax.x - bbox.pMin.x) / dx) + 1;
		dimy = (int)ceil((bbox.pMax.y - bbox.pMin.y) / dy) + 1;

		// allocate data
		typeData = new int[dimx * dimy];
		initData = new CondData2D[dimx * dimy];


		// convert physical coordinates to the grid coordinates
		for (int j = 0; j < num_frames; j++)
			for (int i = 0; i < num_points; i++) 
			{
				points[j][i].x = (points[j][i].x - bbox.pMin.x) / dx;
				points[j][i].y = (points[j][i].y - bbox.pMin.y) / dy;
			}

	}

	void Grid2D::Build(int num_points, Point2D *points, Vec2D* vels)
	{
		// mark all cells as inner 
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				SetType(i, j, IN);
     
		// rasterize lines
        for (int i = 0; i < num_points - 1; i++) 
            RasterLine(points[i], points[i+1], vels[i], vels[i+1], NUM_STEPS, BOUND);
        RasterLine(points[0], points[num_points-1], vels[0], vels[num_points-1], NUM_STEPS, VALVE);

        FloodFill(OUT); 
	}

	void Grid2D::FillTestInitData(Vec2D startVel)
	{
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				switch (GetType(i, j))
				{
				case IN: case OUT: SetData(i, j, CondData2D(NONE, Vec2D(0, 0), 1.0)); break; 
				//case BOUND: SetData(i, j, CondData2D(NOSLIP, Vec2D(0, 0), 1.0)); break; 
				case VALVE: SetData(i, j, CondData2D(FREE, startVel, 1.0)); break; 
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
		fscanf_s(file, "%i", &num_points);

		points = new Point2D*[num_frames];
		velocities = new Vec2D[num_frames];
		Point2D p;
		
		for (int j=0; j<num_frames; j++)
		{
			points[j] = new Point2D[num_points];
			//fscanf_s(file, "");
			for (int i = 0; i < num_points; i++)
				ReadPoint2D(file, points[j][i]);

			ReadPoint2D(file, p);
			velocities[j] = Vec2D(p.x, p.y);
		}
	
		Init();
		return OK;
	}

	Vec2D* Grid2D::ComputeBorderVelocities(int frame)
	{
		Vec2D* res = new Vec2D[num_points];
		int nextframe = (frame == num_frames - 1) ? frame : frame + 1;
		
		for (int i=0; i<num_points; i++)
		{
			res[i].x = (points[nextframe][i].x - points[frame][i].x) * 0.06;
			res[i].y = (points[nextframe][i].y - points[frame][i].y) * 0.06;
		}
		return res;
	}

	void Grid2D::Prepare(int frame)
	{
		if (frame >= num_frames) frame = num_frames-1;
		Vec2D* tempvec = ComputeBorderVelocities(frame);
		Build(num_points, points[frame], tempvec);
		FillTestInitData(velocities[frame]);
		delete[] tempvec;
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