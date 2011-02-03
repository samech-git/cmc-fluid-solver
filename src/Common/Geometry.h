#pragma once

#include <math.h>

#define FTYPE			float
#define INF				1e10
#define COMP_EPS		1e-8
#define BBOX_PADDING	0.02
#define MISSING_VALUE	99999.0f

namespace Common
{
	enum BackendType { CPU, GPU };

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

	enum DirType { X, Y, Z, Z_as_Y };

	struct Vec2D
	{
		FTYPE x, y; 

		Vec2D() : x(0.0), y(0.0) { }
		Vec2D(FTYPE _x, FTYPE _y) : x(_x), y(_y) { }
		Vec2D(Vec2D &vec) : x(vec.x), y(vec.y) { }

		Vec2D operator += (const Vec2D &vec)
		{
			x += vec.x;
			y += vec.y;
			return (*this);
		}

		Vec2D operator * (const FTYPE t)
		{
			return Vec2D(x*t, y*t);
		}

		Vec2D operator + (const Vec2D &vec)
		{
			return Vec2D(x + vec.x, y + vec.y);
		}

		Vec2D operator - (const Vec2D &vec)
		{
			return Vec2D(x - vec.x, y - vec.y);
		}

		FTYPE dot(const Vec2D &vec)
		{
			return (x * vec.x + y * vec.y);
		}

		FTYPE length()
		{
			return sqrt(x*x + y*y);
		}

		void normalize()
		{
			FTYPE t = 1 / length();
			x *= t;
			y *= t;
		}
	};

	struct Vec3D
	{
		FTYPE x, y, z; 

		Vec3D() : x(0.0), y(0.0), z(0.0) { }
		Vec3D(FTYPE _x, FTYPE _y, FTYPE _z) : x(_x), y(_y), z(_z) { }
		Vec3D(Vec3D &vec) : x(vec.x), y(vec.y), z(vec.z) { }

		Vec3D operator += (const Vec3D &vec)
		{
			x += vec.x;
			y += vec.y;
			z += vec.z;
			return (*this);
		}

		Vec3D operator * (const FTYPE t)
		{
			return Vec3D(x*t, y*t, z*t);
		}

		Vec3D operator / (const FTYPE t)
		{
			FTYPE s = 1/t;
			return Vec3D(x*s, y*s, z*s);
		}

		Vec3D operator + (const Vec3D &vec)
		{
			return Vec3D(x + vec.x, y + vec.y, z + vec.z);
		}

		Vec3D operator - (const Vec3D &vec)
		{
			return Vec3D(x - vec.x, y - vec.y, z - vec.z);
		}

		FTYPE length()
		{
			return sqrt(x*x + y*y + z*z);
		}

		bool equal(const Vec3D &vec)
		{
			return ((fabs(x - vec.x) < COMP_EPS) &&
					(fabs(y - vec.y) < COMP_EPS) &&
					(fabs(z - vec.z) < COMP_EPS));
		}

		FTYPE dot(const Vec3D &vec)
		{
			return (x * vec.x + y * vec.y + z * vec.z);
		}

		Vec3D cross(const Vec3D &vec)
		{
			return Vec3D( y * vec.z - z * vec.y, z * vec.x - x * vec.z, x * vec.y - y * vec.x );
		}

		void normalize()
		{
			FTYPE t = 1 / length();
			x *= t;
			y *= t;
			z *= t;
		}
	};

	struct VecTN
	{
		Vec2D tangent, normal;

		VecTN() : tangent(0,0), normal(0,0) { }
		VecTN(Vec2D _x, Vec2D _y) : tangent(_x), normal(_y) { }
		VecTN(VecTN &vec) : tangent(vec.tangent), normal(vec.normal) { }
		VecTN(FTYPE x1, FTYPE y1, FTYPE x2, FTYPE y2) : tangent(x1, y1), normal(x2, y2) { }
	};

	struct Shape2D
	{
		Vec2D* Points;
		Vec2D* Velocities;
		int NumPoints;
		bool Active;

		void Init(int num)
		{
			NumPoints = num;
			Points = new Vec2D[num];
			Velocities = new Vec2D[num];
		}

		void Dispose()
		{
			delete[] Points;
			delete[] Velocities;
		}
	};

	struct Shape3D
	{
		Vec3D* Vertices;
		Vec3D* Velocities;
		int* Indices;

		int NumVertices;
		int NumIndices;
		
		bool Active;

		void Init(Shape3D &shape)
		{
			InitVerts(shape.NumVertices);
			InitInds(shape.NumIndices);
		}

		void InitVerts(int num)
		{
			NumVertices = num;
			Vertices = new Vec3D[num];
			Velocities = new Vec3D[num];
		}

		void InitInds(int num)
		{
			NumIndices = num;
			Indices = new int[num*3];
		}

		Shape3D() : Vertices(0), Velocities(0), Indices(0) { }
	
		~Shape3D()
		{
			if (Vertices) delete[] Vertices;
			if (Velocities) delete[] Velocities;
			if (Indices) delete[] Indices;
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
			return ((fabs(MinX - dest.MinX) < COMP_EPS) && 
					(fabs(MinY - dest.MinY) < COMP_EPS) && 
					(fabs(Dx - dest.Dx) < COMP_EPS) && 
					(fabs(Dy - dest.Dy) < COMP_EPS) && 
					Nx == dest.Nx && 
					Ny == dest.Ny);
		}

		Vec2D GetVelocity(double x, double y)
		{
			double tx = (x - MinX) / Dx;
			double ty = (y - MinY) / Dy;
			
			if (tx<0 || ty<0 || tx>=Nx-1 || ty>=Ny-1 || x==0)
				return Vec2D(0, 0);

			int itx = (int)tx;
			int ity = (int)ty;

			int t = itx + ity*Nx;
			
			FTYPE rx = Data[t].x;
			FTYPE ry = Data[t].y;

			return Vec2D(rx, ry);
		}
	};

	struct Field3D
	{
		double MinX, MinY, MinZ;
		double Dx, Dy, Dz;
		int Nx, Ny, Nz;
		Vec3D* Data;

		void Init()
		{
			MinX = MinY = MinZ = 0.0;
			Dx = Dy = Dz = 0.0;
			Nx = Ny = Nz = 0;
			Data = 0;
		}

		void Init(double minx, double miny, double minz, double dx, double dy, double dz, int nx, int ny, int nz)
		{
			MinX = minx; MinY = miny; MinZ = minz;
			Dx = dx; Dy = dy; Dz = dz;
			Nx = nx; Ny = ny; Nz = nz;
			Data = new Vec3D[nx * ny * nz];
		}

		void Init(Field3D src)
		{
			Init(src.MinX, src.MinY, src.MinZ, src.Dx, src.Dy, src.Dz, src.Nx, src.Ny, src.Nz);
		}

		void Dispose()
		{
			if (Data != 0) delete[] Data;
		}

		bool Correlate(Field3D dest)
		{
			return ((fabs(MinX - dest.MinX) < COMP_EPS) && 
					(fabs(MinY - dest.MinY) < COMP_EPS) && 
					(fabs(MinZ - dest.MinZ) < COMP_EPS) && 
					(fabs(Dx - dest.Dx) < COMP_EPS) && 
					(fabs(Dy - dest.Dy) < COMP_EPS) && 
					(fabs(Dz - dest.Dz) < COMP_EPS) && 
					Nx == dest.Nx && 
					Ny == dest.Ny &&
					Nz == dest.Nz);
		}

		Vec3D GetVelocity(double x, double y, double z)
		{
			double tx = (x - MinX) / Dx;
			double ty = (y - MinY) / Dy;
			double tz = (z - MinZ) / Dz;
			
			if (tx<0 || ty<0 || tz<0 || tx>=Nx-1 || ty>=Ny-1 || tz>=Nz-1)
				return Vec3D(0, 0, 0);

			int itx = (int)tx;
			int ity = (int)ty;
			int itz = (int)tz;

			int t = itx + ity*Nx + itz*Nx*Ny;
			
			double rx = Data[t].x;
			double ry = Data[t].y;
			double rz = Data[t].z;

			return Vec3D((FTYPE)rx, (FTYPE)ry, (FTYPE)rz);
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

	struct FrameInfo3D
	{
		Shape3D* Shapes;
		int NumShapes;
		double Duration;
		
		Field3D Field;

		void Init(int num)
		{
			NumShapes = num;
			Shapes = new Shape3D[num];
			Field.Init();
		}

		~FrameInfo3D()
		{
			delete[] Shapes;
			Field.Dispose();
		}
	};

	struct DepthInfo3D 
	{
		int dimx, dimy;
		float *depth;

		DepthInfo3D(int nx, int ny)
		{
			dimx = nx; dimy = ny;
			depth = new float[nx * ny];
		};

		~DepthInfo3D()
		{
			delete [] depth;
		}
	};

	struct BBox2D
	{
		Vec2D pMin, pMax;

		BBox2D() { Clear(); }
		
		void AddPoint(Vec2D p)
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

			FTYPE wx = pMax.x - pMin.x;
			FTYPE wy = pMax.y - pMin.y;
			
			pMin.x -= wx * (FTYPE)BBOX_PADDING;
			pMin.y -= wy * (FTYPE)BBOX_PADDING;
			
			pMax.x += wx * (FTYPE)BBOX_PADDING;
			pMax.y += wy * (FTYPE)BBOX_PADDING;
		}

		void Clear()
		{
			pMin.x = pMin.y = INF; 
			pMax.x = pMax.y = -INF;
		}
	};

	struct BBox3D
	{
		Vec3D pMin, pMax;

		BBox3D() { Clear(); }
		
		void AddPoint(Vec3D p)
		{
			if (p.x < pMin.x) pMin.x = p.x;
			if (p.y < pMin.y) pMin.y = p.y;
			if (p.z < pMin.z) pMin.z = p.z;
			
			if (p.x > pMax.x) pMax.x = p.x;
			if (p.y > pMax.y) pMax.y = p.y;
			if (p.z > pMax.z) pMax.z = p.z;
		}

		void Build(int num_frames, FrameInfo3D* frames)
		{
			Clear();
			for (int j = 0; j < num_frames; j++)
				for (int i = 0; i < frames[j].NumShapes; i++)
					for (int k = 0; k < frames[j].Shapes[i].NumVertices; k++)
						AddPoint(frames[j].Shapes[i].Vertices[k]);

			FTYPE wx = pMax.x - pMin.x;
			FTYPE wy = pMax.y - pMin.y;
			FTYPE wz = pMax.z - pMin.z;
			
			pMin.x -= wx * (FTYPE)BBOX_PADDING;
			pMin.y -= wy * (FTYPE)BBOX_PADDING;
			pMin.z -= wz * (FTYPE)BBOX_PADDING;
			
			pMax.x += wx * (FTYPE)BBOX_PADDING;
			pMax.y += wy * (FTYPE)BBOX_PADDING;
			pMax.z += wz * (FTYPE)BBOX_PADDING;
		}

		void Clear()
		{
			pMin.x = pMin.y = pMin.z = INF; 
			pMax.x = pMax.y = pMax.z = -INF;
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