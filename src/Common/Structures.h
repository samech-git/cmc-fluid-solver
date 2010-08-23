#pragma once

#define FTYPE			float
#define INF				1e10

namespace Common
{
	struct Vec2D
	{
		double x, y; 

		Vec2D() : x(0.0), y(0.0) { }
		Vec2D(double _x, double _y) : x(_x), y(_y) { }
		Vec2D(Vec2D &vec) : x(vec.x), y(vec.y) { }
	};

	struct Vec3D
	{
		FTYPE x, y, z; 

		Vec3D() : x(0.0), y(0.0), z(0.0) { }
		Vec3D(FTYPE _x, FTYPE _y, FTYPE _z) : x(_x), y(_y), z(_z) { }
		Vec3D(Vec3D &vec) : x(vec.x), y(vec.y), z(vec.z) { }
	};

	struct VecTN
	{
		Vec2D tangent, normal;

		VecTN() : tangent(0,0), normal(0,0) { }
		VecTN(Vec2D _x, Vec2D _y) : tangent(_x), normal(_y) { }
		VecTN(VecTN &vec) : tangent(vec.tangent), normal(vec.normal) { }
		VecTN(double x1, double y1, double x2, double y2) : tangent(x1, y1), normal(x2, y2) { }
	};

	struct Point2D 
	{ 
		double x, y; 

		Point2D() : x(0.0), y(0.0) { }
		Point2D(double _x, double _y) : x(_x), y(_y) { }
	};

	struct Shape2D
	{
		Point2D* Points;
		Vec2D* Velocities;
		int NumPoints;
		bool Active;

		void Init(int num)
		{
			NumPoints = num;
			Points = new Point2D[num];
			Velocities = new Vec2D[num];
		}

		void Dispose()
		{
			delete[] Points;
			delete[] Velocities;
		}
	};


	struct Field2D
	{
		double MinX;
		double MinY;
		double Dx;
		double Dy;
		int Nx;
		int Ny;
		Vec2D* Data;

		void Init()
		{
			MinX = 0;
			MinY = 0;
			Dx = 0;
			Dy = 0;
			Nx = 0;
			Ny = 0;
			Data = 0;
		}

		void Init(double minx, double miny, double dx, double dy, int nx, int ny)
		{
			MinX = minx;
			MinY = miny;
			Dx = dx;
			Dy = dy;
			Nx = nx;
			Ny = ny;
			Data = new Vec2D[nx*ny];
		}

		void Init(Field2D src)
		{
			Init(src.MinX, src.MinY, src.Dx, src.Dy, src.Nx, src.Ny);
		}

		void Dispose()
		{
			if (Data != 0) delete[] Data;
		}

		bool Correlate(Field2D dest)
		{
			return (MinX == dest.MinX && MinY == dest.MinY && Dx == dest.Dx && Dy == dest.Dy && Nx == dest.Nx && Ny == dest.Ny);
		}

		Vec2D GetVelocity(double x, double y)
		{
			double tx = (x - MinX) / Dx;
			double ty = (y - MinY) / Dy;
			
			if (tx<0 || ty<0 || tx>=Nx-1 || ty>=Ny-1)
				return Vec2D(0, 0);

			int itx = (int)tx;
			int ity = (int)ty;
			/*double dtx = tx - itx;
			double dty = ty - ity;*/

			/*double f1 = (1 - dtx) * (1 - dty);
			double f2 = dtx * (1 - dty);
			double f3 = (1 - dtx) * dty;
			double f4 = dtx * dty;*/

			int t = itx + ity*Nx;
			
			double rx = Data[t].x;
			double ry = Data[t].y;

			return Vec2D(rx, ry);
		}
	};

	struct FrameInfo2D
	{
		Shape2D* Shapes;
		int NumShapes;
		double Duration;
		
		Field2D Field;

		void Init(int num)
		{
			NumShapes = num;
			Shapes = new Shape2D[num];
			Field.Init();
		}
		void Dispose()
		{
			for (int i=0; i<NumShapes; i++)
				Shapes[i].Dispose();
			delete[] Shapes;
			Field.Dispose();
		}
	};


	struct BBox2D
	{
		Point2D pMin, pMax;

		BBox2D() { Clear(); }
		
		void AddPoint(Point2D p)
		{
			if (p.x < pMin.x) pMin.x = p.x;
			if (p.y < pMin.y) pMin.y = p.y;
			if (p.x > pMax.x) pMax.x = p.x;
			if (p.y > pMax.y) pMax.y = p.y;
		}

		void Build(int num_frames, FrameInfo2D* frames)
		{
			Clear();
			for (int j = 0; j < num_frames; j++)
				for (int i = 0; i < frames[j].NumShapes; i++)
					for (int k = 0; k < frames[j].Shapes[i].NumPoints; k++)
						AddPoint(frames[j].Shapes[i].Points[k]);

			double wx = pMax.x - pMin.x;
			double wy = pMax.y - pMin.y;
			
			pMin.x -= wx * 0.02;
			pMin.y -= wy * 0.02;
			
			pMax.x += wx * 0.02;
			pMax.y += wy * 0.02;
		}

		void Clear()
		{
			pMin.x = pMin.y = INF; 
			pMax.x = pMax.y = -INF;
		}
	};

	struct FluidParams
	{
		FTYPE v_T, v_vis;
		FTYPE t_vis, t_phi;

		FluidParams() { }
		
		FluidParams(double Re, double Pr, double lambda)  
		{
			v_T = 1.0;
			v_vis = (FTYPE)(1.0 / Re);

			t_vis = (FTYPE)(1.0 / (Re * Pr));
			t_phi = (FTYPE)((lambda - 1) / (lambda * Re));
		}

		FluidParams(double vis, double rho, double R, double k, double cv)
		{
			v_T = (FTYPE)R;
			v_vis = (FTYPE)(vis / rho);

			t_vis = (FTYPE)(k / (rho * cv));
			t_phi = (FTYPE)(vis / (rho * cv));
		}
	};

	static int AlignBy32(int num)
	{
		if( (num & 31) == 0 ) return num;
			else return (((num >> 5) + 1) << 5);
	}
}