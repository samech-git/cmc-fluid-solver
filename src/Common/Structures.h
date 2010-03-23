#pragma once

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
		double x, y, z; 

		Vec3D() : x(0.0), y(0.0), z(0.0) { }
		Vec3D(double _x, double _y, double _z) : x(_x), y(_y), z(_z) { }
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

	struct FrameInfo2D
	{
		Shape2D* Shapes;
		int NumShapes;
		double Duration;

		void Init(int num)
		{
			NumShapes = num;
			Shapes = new Shape2D[num];
		}
		void Dispose()
		{
			for (int i=0; i<NumShapes; i++)
				Shapes[i].Dispose();
			delete[] Shapes;
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
		double v_T, v_vis;
		double t_vis, t_phi;

		FluidParams() { }
		
		FluidParams(double Re, double Pr, double lambda)  
		{
			v_T = 1.0;
			v_vis = 1.0 / Re;

			t_vis = 1.0 / (Re * Pr);
			t_phi = (lambda - 1) / (lambda * Re);
		}

		FluidParams(double vis, double rho, double R, double k, double cv)
		{
			v_T = R;
			v_vis = vis / rho;

			t_vis = k / (rho * cv);
			t_phi = vis / (rho * cv);
		}
	};
}