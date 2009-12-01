#include "Grid2D.h"

namespace FluidSolver
{
	Grid2D::Grid2D(double _dx, double _dy) : dx(_dx), dy(_dy), data(NULL)
	{
	}

	Grid2D::~Grid2D()
	{
		if (data != NULL) delete [] data;
	}

	inline CellType Grid2D::GetCell(int x, int y)
	{
		return (CellType)data[x * (dimy + 1) + y];
	}

	inline void Grid2D::SetCell(int x, int y, CellType t)
	{
		data[x * (dimy + 1) + y] = t;
	}

	void ReadPoint2D(FILE *file, Point2D &p)
	{
		string str = "";
		char c;	

		// read line
		fscanf_s(file, "%c", &c);
		while (c == '\n') fscanf_s(file, "%c", &c);
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

	void Grid2D::BuildBBox(int num_points, Point2D *points)
	{
		bbox.Clear();
		for (int i = 0; i < num_points; i++)
			bbox.AddPoint(points[i]);

		bbox.pMin.x -= BBOX_OFFSET;
		bbox.pMin.y -= BBOX_OFFSET;
		bbox.pMax.x += BBOX_OFFSET;
		bbox.pMax.y += BBOX_OFFSET;

#ifdef _DEBUG
		printf("bbox x range: [%f, %f]\n", bbox.pMin.x, bbox.pMax.x);
		printf("bbox y range: [%f, %f]\n", bbox.pMin.y, bbox.pMax.y);
#endif
	}

	void Grid2D::RasterLine(Point2D p1, Point2D p2, int steps, CellType color)
    {
        Point2D dv((p2.x - p1.x) / steps, (p2.y - p1.y) / steps);		// divide a segment into parts
		Point2D startp = p1;

		// go through the whole segment
        for (int i = 0; i <= steps; i++) 
        {
            int x = (int)startp.x;
            int y = (int)startp.y;

			SetCell(x, y, color);
			SetCell(x+1, y, color);
			SetCell(x, y+1, color);
			SetCell(x+1, y+1, color);

			startp.x += dv.x;
            startp.y += dv.y;
        }
    }

	void Grid2D::FloodFill(CellType color)
    {
		// go through all the rows
        for (int j = 0; j <= dimy; j++) 
        {
			// left to right
            int i = 0;
            while (i <= dimx && GetCell(i, j) == 1) 
            {
                SetCell(i, j, color);
                i++;
            }

			// right to left
            i = dimx;
            while (GetCell(i, j) == 1 && i >= 0) 
            {
                SetCell(i, j, color);
                i--;
            }
        }
    }

	void Grid2D::Build(int num_points, Point2D *points)
	{
		BuildBBox(num_points, points);
	
		dimx = (int)ceil((bbox.pMax.x - bbox.pMin.x) / dx);
		dimy = (int)ceil((bbox.pMax.y - bbox.pMin.y) / dy);

		// allocate data
		data = new int[(dimx + 1) * (dimy + 1)];

		// mark all cells as inner 
		for (int i = 0; i <= dimx; i++)
			for (int j = 0; j <= dimy; j++)
				SetCell(i, j, IN);

		// convert physical coordinates to the grid coordinates
		for (int i = 0; i < num_points; i++) 
		{
            points[i].x = (points[i].x - bbox.pMin.x) / dx;
            points[i].y = (points[i].y - bbox.pMin.y) / dy;
        }
     
		// rasterize lines
        for (int i = 0; i < num_points - 1; i++) 
            RasterLine(points[i], points[i+1], NUM_STEPS, BOUND);
        RasterLine(points[0], points[num_points-1], NUM_STEPS, VALVE);

        FloodFill(OUT); 
	}

	int Grid2D::LoadFromFile(char *filename)
	{
		FILE *file = NULL;
		if (fopen_s(&file, filename, "r") != 0)
		{
			printf("Error: cannot open file \"%s\" \n", filename);
			return ERR;
		}

		int num_frames;
		fscanf_s(file, "%i", &num_frames);	// currently not used
		
		// read points
		int num_points;
		fscanf_s(file, "%i", &num_points);
		Point2D *points = new Point2D[num_points];
		for (int i = 0; i < num_points; i++)
			ReadPoint2D(file, points[i]);

		Build(num_points, points);

		delete [] points;

		return OK;
	}

	void Grid2D::TestPrint()
	{
		printf("%i %i\n", dimx, dimy);
		for (int i = 0; i < dimx; i++)
		{
			for (int j = 0; j < dimy; j++)
			{
				CellType t = GetCell(i, j);
				switch (t)
				{
				case IN: printf("."); break;
				case OUT: printf(" "); break;
				case BOUND: printf("#"); break;
				case VALVE: printf("+"); break;
				}
			}
			printf("\n");
		}
	}
}