/*
 *  Copyright 2010-2011 Nikolai Sakharnykh
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include "Grid3D.h"

#ifdef linux
#include <cmath> //for abs functions
#endif

namespace FluidSolver3D
{
	Grid3D::Grid3D(double _dx, double _dy, double _dz, double _depth, double _baseT, BackendType _backend, bool useNetCDF, SplitType _split_type) : 
		dx(_dx), dy(_dy), dz(_dz), depth(_depth), baseT(_baseT), nodes(NULL), d_types(NULL), d_typesT(NULL), backend(_backend), use3Dshape(false), frames(NULL), depthInfo(NULL), split_type(_split_type)
	{
		grid2D = new FluidSolver2D::Grid2D(dx, dy, baseT, true, 0.0);
	}

	Grid3D::Grid3D(double _dx, double _dy, double _dz, double _baseT, BackendType _backend, bool _useNetCDF, SplitType _split_type) : 
		dx(_dx), dy(_dy), dz(_dz), baseT(_baseT), nodes(NULL), d_types(NULL), d_typesT(NULL), backend(_backend), use3Dshape(true), useNetCDF(_useNetCDF), frames(NULL), depthInfo(NULL), split_type(_split_type)
	{
		grid2D = NULL;
	}

	Grid3D::~Grid3D()
	{
		if (nodes != NULL) delete [] nodes;
		if (d_typesT != NULL) multiDevFree<NodeType>(d_typesT);
		if (d_types != NULL) multiDevFree<NodeType>(d_types);
		if (grid2D != NULL) delete grid2D;
		if (frames != NULL) delete [] frames;
		if (depthInfo != NULL) delete depthInfo;
	}

	void Grid3D::GenerateListSegments(int &numSeg, Segment3D *h_list, int dim1, int dim2, int dim3, DirType dir, int nblockZ = 1)
	{
		numSeg = 0;
		int _nblockZ = 1;
		switch (dir)
		{
		case X:
		case Y:
			_nblockZ = nblockZ;
			break;
		}
		int block_start = 0;
		int block_end = dim3 / _nblockZ; // will sort by blocks along Z

		while (block_end <= dim3)
		{
			for (int i = 0; i < dim2; i++)
				for (int j = block_start; j < block_end; j++)
				{
					Segment3D seg, new_seg;
					int state = 0, incx, incy, incz;
					switch (dir)
					{
					case X:	
						seg.posx = 0; seg.posy = i; seg.posz = j; 
						incx = 1; incy = 0; incz = 0;
						break;
					case Y: 
						seg.posx = i; seg.posy = 0; seg.posz = j; 
						incx = 0; incy = 1; incz = 0;
						break;
					case Z: 
						seg.posx = i; seg.posy = j; seg.posz = 0; 
						incx = 0; incy = 0; incz = 1; 
						break;
					}
					seg.dir = (DirType)dir;

					while ((seg.posx + incx < dimx) && (seg.posy + incy < dimy) && (seg.posz + incz < dimz))
					{
						if (GetType(seg.posx + incx, seg.posy + incy, seg.posz + incz) == NODE_IN)
						{
							if (state == 0) 
								new_seg = seg;
							state = 1;
						}
						else
						{
							if (state == 1)
							{
								new_seg.endx = seg.posx + incx;
								new_seg.endy = seg.posy + incy;
								new_seg.endz = seg.posz + incz;

								new_seg.size = (new_seg.endx - new_seg.posx) + (new_seg.endy - new_seg.posy) + (new_seg.endz - new_seg.posz) + 1;

								new_seg.skipX = false;
								new_seg.type = BOUND;
								h_list[numSeg] = new_seg;
								numSeg++;
								state = 0;
							}
						}
					
						seg.posx += incx;
						seg.posy += incy;
						seg.posz += incz;
					}
				}
				block_start = block_end;
				block_end += (block_end + dim3 / _nblockZ == (dim3 / _nblockZ)* _nblockZ)? dim3 / _nblockZ + dim3 % _nblockZ: dim3 / _nblockZ;
		}
	}

	void Grid3D::GenerateGridBoundaries(NodesBoundary3D *node_list, int numSeg, Segment3D *h_list, bool transposed)
	{
			for (int iseg = 0; iseg < numSeg; iseg++)
			{
				node_list[iseg].first = nodes[ h_list[iseg].posx * dimy * dimz + h_list[iseg].posy * dimz + h_list[iseg].posz ];
				node_list[iseg].last = nodes[ h_list[iseg].endx * dimy * dimz + h_list[iseg].endy * dimz + h_list[iseg].endz ];
				//if (transposed)
				//{
				//	node_list[iseg].first = nodes[ h_list[iseg].posx * dimy * dimz + h_list[iseg].posy * dimz + h_list[iseg].posz ];
				//	node_list[iseg].last = nodes[ h_list[iseg].endx * dimy * dimz + h_list[iseg].endy * dimz + h_list[iseg].endz ];
				//}
				//else
				//{
				//	node_list[iseg].first = nodes[ h_list[iseg].posx * dimy * dimz + h_list[iseg].posz * dimy + h_list[iseg].posy ];
				//	node_list[iseg].last = nodes[ h_list[iseg].endx * dimy * dimz + h_list[iseg].endz * dimy + h_list[iseg].endy ];
				//}
			}
	}

	void Grid3D::SplitSegments_X(int *splitting)
	{
		if (backend == CPU) // not supported in MPI 
		{
			splitting[0] = dimx;
			return;
		}
		int i, num_seg_X, num_seg_Y, num_seg_Z;
		double loadPerGPU;
		double *acu_sum = new double[dimx];

		Segment3D *h_list_X = NULL, *h_list_Y = NULL, *h_list_Z = NULL;

		PARAplan *pplan = PARAplan::Instance();
		int nGPUs = pplan->gpuTotal();

		for (i = 0; i < dimx; i++)
			acu_sum[i] = 0.0;

		if (split_type == EVEN_SEGMENTS)
		{
			h_list_X = new Segment3D[dimy * dimz * MAX_SEGS_PER_ROW];
			h_list_Y = new Segment3D[dimx * dimz * MAX_SEGS_PER_ROW];
			h_list_Z = new Segment3D[dimx * dimy * MAX_SEGS_PER_ROW];

			GenerateListSegments(num_seg_X, h_list_X, dimx, dimy, dimz, X);	
			GenerateListSegments(num_seg_Y, h_list_Y, dimy, dimx, dimz, Y);	
			GenerateListSegments(num_seg_Z, h_list_Z, dimz, dimx, dimy, Z);	

			for (i = 0; i < num_seg_Y; i++)
				acu_sum[h_list_Y[i].posx] += 1.0;
			for (i = 0; i < num_seg_Z; i++)
				acu_sum[h_list_Z[i].posx] += 1.0;
		
			for (int ii = 0; ii < num_seg_X; ii++)
				for (i = h_list_X[ii].posx; i <= h_list_X[ii].endx; i++)
					acu_sum[i] += 1.0/h_list_X[ii].size;

			loadPerGPU = double(num_seg_Y + num_seg_Z + num_seg_X) / nGPUs;
		}

		if (split_type == EVEN_VOLUME)
		{
			double indsidePoints = 0.;
			for (i = 0; i < dimx; i++)
				for (int j = 0; j < dimy; j++)
					for (int k = 0; k < dimz; k++)
						if (GetType(i,j,k) == NODE_IN)
						{
							acu_sum[i] += 1.0;
							indsidePoints += 1.0;
						}
			loadPerGPU = indsidePoints / nGPUs;
		}

		// save sum to file:
		
		if (pplan->rank() == 0)
		{
			FILE *file = NULL;
			fopen_s(&file, "segs_per_X.txt", "w");
			for (i = 0; i < dimx; i++) 
				fprintf(file, "%i    %f\n",i , acu_sum[i]);
			fclose(file);
		}
		
		double sum = acu_sum[0];
		int gpunode = 0;
		int i_old = 0 ;
		for (i = 1; i < dimx; i++)
		{
			if (sum + acu_sum[i] > loadPerGPU)
			{
				splitting[gpunode] = i - 1 - i_old + 1;
				gpunode++;
				i_old = i;
				if (gpunode >= nGPUs - 1) break;				
				sum = 0.0;
			}
			sum += acu_sum[i];	
		}
		splitting[nGPUs - 1] = (nGPUs == 1)? dimx: dimx - i_old;

		if (h_list_Z != NULL) delete [] h_list_Z;
		if (h_list_Y != NULL) delete [] h_list_Y;
		if (h_list_X != NULL) delete [] h_list_X;
		delete [] acu_sum;
	}

	BBox3D Grid3D::GetBBox()
	{
		return bbox;
	}

	NodeType Grid3D::GetType(int i, int j, int k)
	{
		return nodes[i * dimy * dimz + j * dimz + k].type;
	}

	NodeType **Grid3D::GetTypesGPU(bool transposed)
	{
		return (transposed ? d_typesT : d_types);
	}

	Node *Grid3D::GetNodesCPU()
	{
		return nodes;
	}

	BCtype Grid3D::GetBC_vel(int i, int j, int k)
	{
		return nodes[i * dimy * dimz + j * dimz + k].bc_vel;
	}

	BCtype Grid3D::GetBC_temp(int i, int j, int k)
	{
		return nodes[i * dimy * dimz + j * dimz + k].bc_temp;
	}

	Vec3D Grid3D::GetVel(int i, int j, int k)
	{
		return nodes[i * dimy * dimz + j * dimz + k].v;
	}

	FTYPE Grid3D::GetT(int i, int j, int k)
	{
		return nodes[i * dimy * dimz + j * dimz + k].T;
	}

	void Grid3D::SetType(int i, int j, int k, NodeType _type)
	{
		nodes[i * dimy * dimz + j * dimz + k].type = _type;
	}

	void Grid3D::SetData(int i, int j, int k, BCtype _bc_vel, BCtype _bc_T, const Vec3D &_vel, FTYPE _T)
	{
		int index = i * dimy * dimz + j * dimz + k;
		nodes[index].bc_vel = _bc_vel;
		nodes[index].bc_temp = _bc_T;
		nodes[index].v = _vel;
		nodes[index].T = _T;
	}

	double Grid3D::GetFrameTime()
	{
		int frames = grid2D->GetFramesNum();
		double length = grid2D->GetCycleLenght();
		return (length / frames);
	}

	int Grid3D::GetFramesNum()
	{
		return num_frames;
	}

	double Grid3D::GetCycleLength()
	{
		if( use3Dshape )
			return frame_time;
		else
			return grid2D->GetCycleLenght();
	}

	void Grid3D::SetFrameTime(double time)
	{
		frame_time = time;
	}

	void Grid3D::SetBoundParams(const Vec3D &vec, const double &temp)
	{
		bcInVel = vec;
		bcInT = temp;
	}

	int Grid3D::GetFrame(double time)
	{
		if( use3Dshape )
			return 0;
		else
			return grid2D->GetFrame(time);
	}
	
	float Grid3D::GetLayerTime(double time)
	{
		if( use3Dshape )
			return (float)frame_time;
		else
			return grid2D->GetLayerTime(time);
	}

	FluidSolver2D::Grid2D *Grid3D::GetGrid2D()
	{
		return grid2D;
	}

	DepthInfo3D *Grid3D::GetDepthInfo()
	{
		return depthInfo;
	}

	void Grid3D::SetNodeVel(int i, int j, int k, Vec3D new_v)
	{
		nodes[i * dimy * dimz + j * dimz + k].v = new_v;
	}

	void Grid3D::Init(bool align)
	{
		dimx = ((int)ceil((bbox.pMax.x - bbox.pMin.x) / dx) + 1);
		dimy = ((int)ceil((bbox.pMax.y - bbox.pMin.y) / dy) + 1);
		dimz = ((int)ceil((bbox.pMax.z - bbox.pMin.z) / dz) + 1);

		if( align ) { dimx = AlignBy32(dimx); dimy = AlignBy32(dimy); dimz = AlignBy32(dimz); } 

		// allocate data
		int size = dimx * dimy * dimz;
		nodes = new Node[size];
		
		for (int i=0; i<dimx * dimy * dimz; i++)
		{
			nodes[i].type = NODE_OUT;
			nodes[i].v = Vec3D(0.0, 0.0, 0.0);
			nodes[i].T = 0;
		}
	}

	bool Grid3D::Load3DShape(char *filename, bool align)
	{
		// load 3D poly model
		FILE *file = NULL;
		if (fopen_s(&file, filename, "r") != 0)
		{
			printf("Error: cannot open file \"%s\" \n", filename);
			return false;
		}

		fscanf_s(file, "%i", &num_frames);	
		Vec3D p;
		int temp;
		frames = new FrameInfo3D[num_frames];
	
		for (int j=0; j<num_frames; j++)
		{
			// use only 1 shape for now
			frames[j].Init(1);

			for (int i = 0; i<frames[j].NumShapes; i++)
			{
				fscanf_s(file, "%i", &temp);
				frames[j].Shapes[i].InitVerts(temp);
				for (int k=0; k<frames[j].Shapes[i].NumVertices; k++)
				{
					ReadPoint3D(file, p);
					frames[j].Shapes[i].Vertices[k] = p * GRID_SCALE_FACTOR;
					
					ReadPoint3D(file, p);
					frames[j].Shapes[i].Velocities[k] = p;
				}

				fscanf_s(file, "%i", &temp);
				frames[j].Shapes[i].InitInds(temp);
				for (int k=0; k<frames[j].Shapes[i].NumIndices*3; k++)
				{
					fscanf_s(file, "%i", &temp);
					frames[j].Shapes[i].Indices[k] = temp;
				}

				frames[j].Shapes[i].Active = false;
				frames[j].Duration = 1.0 / 75;			// 75 fps
			}
		}
		fclose(file);

		bbox.Build(num_frames, frames);
		
		Init(align);

		// convert physical coordinates to the grid coordinates
		for (int j = 0; j < num_frames; j++)
			for (int i = 0; i < frames[j].NumShapes; i++)
				for (int k = 0; k < frames[j].Shapes[i].NumVertices; k++)
				{
					frames[j].Shapes[i].Vertices[k].x = (frames[j].Shapes[i].Vertices[k].x - bbox.pMin.x) / (FTYPE)dx;
					frames[j].Shapes[i].Vertices[k].y = (frames[j].Shapes[i].Vertices[k].y - bbox.pMin.y) / (FTYPE)dy;
					frames[j].Shapes[i].Vertices[k].z = (frames[j].Shapes[i].Vertices[k].z - bbox.pMin.z) / (FTYPE)dz;
				}
		
		return true;
	}

	bool Grid3D::LoadNetCDF(char *filename, bool align)
	{
		int status;

		// open netcdf
		int ncid;
		status = nc_open(filename, NC_NOWRITE, &ncid);

		// read dimensions
		int latid, lonid;
		nc_inq_dimid(ncid, "_lat_subset", &latid);
		nc_inq_dimid(ncid, "_lon_subset", &lonid);

		size_t nx, ny;
		nc_inq_dimlen(ncid, latid, &nx);
		nc_inq_dimlen(ncid, lonid, &ny);

		double *lats = new double[nx];
		nc_get_var(ncid, latid, lats);

		double *lons = new double[ny];
		nc_get_var(ncid, lonid, lons);
	
		// read depths
		int varid;
		nc_inq_varid(ncid, "z", &varid);

		depthInfo = new DepthInfo3D(nx, ny);
		nc_get_var(ncid, varid, depthInfo->depth);

		// build bbox
		bbox.AddPoint(Vec3D((float)lats[0], (float)lons[0], 0.0f));
		bbox.AddPoint(Vec3D((float)lats[nx-1], (float)lons[ny-1], 0.0f));
		for (int j = 0; j < (int)ny; j++)
			for (int i = 0; i < (int)nx; i++)
			{
				float z = depthInfo->depth[i + j * nx];
				if( z < bbox.pMin.z ) bbox.pMin.z = z;
			}		
		bbox.pMin.z -= (float)dz;

		delete [] lats;
		delete [] lons;

		Init(align);

		num_frames = 1;

		return true;
	}

	bool Grid3D::LoadFromFile(char *filename, bool align)
	{
		if (use3Dshape)
		{
			if (useNetCDF) 
				return LoadNetCDF(filename, align);
			else 
				return Load3DShape(filename, align);
		}
		else
		{
			// load 2D shape, extend in depth
			if (grid2D->LoadFromFile(filename, "", align))
			{
				dimx = grid2D->dimx;
				dimy = grid2D->dimy;
				active_dimz = (int)ceil(depth / dz) + 1;
				dimz = align ? AlignBy32(active_dimz) : active_dimz;
				nodes = new Node[dimx * dimy * dimz];
				num_frames = grid2D->GetFramesNum();
				return true;
			}
			else
				return false;
		}
	}

	void Grid3D::Prepare_CPU(double time)
	{
		if (use3Dshape) 
		{
			if (useNetCDF) Prepare3D_NetCDF(time);
				else Prepare3D_Shape(time);
		}
		else
			Prepare2D(time);
	}

	void Grid3D::Prepare_GPU()
	{
		PARAplan* pplan = PARAplan::Instance();
		int dimxOffset = pplan->getOffset1D();
		int dimxNode = pplan->getLength1D();

		NodeType *nodeTypes = new NodeType[dimxNode * dimy * dimz];
		NodeType *nodeTypesT = new NodeType[dimxNode * dimy * dimz];

		for (int i = dimxOffset; i < dimxOffset + dimxNode; i++)
			for (int j = 0; j < dimy; j++)
				for (int k = 0; k < dimz; k++)
				{
					int id = (i) * dimy * dimz + j * dimz + k;
					int idT = (i) * dimy * dimz + k * dimy + j;
					nodeTypes[id - dimxOffset*dimy*dimz] = nodes[id].type;
					nodeTypesT[id - dimxOffset*dimy*dimz] = nodes[idT].type;
				}
#if (TRANSPOSE_OPT == 1)		
		// currently implemented on CPU but it's possible to do a transpose on GPU
		//for (int i = dimxOffset; i < dimxOffset + dimxNode; i++)
		//	for (int j = 0; j < dimy; j++)
		//		for (int k = j+1; k < dimz; k++)
		//		{
		//			int id = (i) * dimy * dimz + j * dimz + k;
		//			int idT = (i) * dimy * dimz + k * dimy + j;
		//			Node tmp = nodes[idT];
		//			nodes[idT] = nodes[id];
		//			nodes[id] = tmp;
		//		}
		//multiDevMemcpy<Node>(d_nodesT, nodes + dimxOffset*dimy*dimz, dimxNode * dimy * dimz);
		//multiDevMemcpy<Node>(nodes + dimxOffset*dimy*dimz, d_nodes, dimxNode * dimy * dimz);

		multiDevMemcpy<NodeType>(d_types, nodeTypes, dimxNode * dimy * dimz);
		multiDevMemcpy<NodeType>(d_typesT, nodeTypesT, dimxNode * dimy * dimz);
#endif

		delete [] nodeTypesT;
		delete [] nodeTypes;
	}

	void Grid3D::Split()
	{
		PARAplan *pplan = PARAplan::Instance();
		switch (split_type)
		{
		case EVEN_X:
			pplan->splitEven1D(dimx);
			break;
		default:
			int *splitting = new int[dimx];
			SplitSegments_X(splitting);
			pplan->split1D(splitting);
			delete [] splitting;
		}
	}

	void Grid3D::Init_GPU()
	{
		if (backend != GPU) return;
		
		if (d_typesT != NULL) multiDevFree<NodeType>(d_typesT);
		if (d_types != NULL) multiDevFree<NodeType>(d_types);
		PARAplan *pplan = PARAplan::Instance();
		int size = pplan->getLength1D() * dimy * dimz;
#if (TRANSPOSE_OPT == 1)
#endif
		multiDevAlloc<NodeType> (d_types, size); 
		multiDevAlloc<NodeType> (d_typesT, size); 
		Prepare_GPU();
	}

	void Grid3D::Prepare(double time)
	{	
		Prepare_CPU(time);
		// copy to GPU as well
		if( backend == GPU )
		{
			Prepare_GPU();
		}
	}

	void Grid3D::Prepare2D(double time)
	{
		grid2D->Prepare(time);

		memset(nodes, 0, sizeof(Node) * dimx * dimy * dimz);
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
			{
				int ind;
				if (grid2D->GetType(i, j) == NODE_OUT)
				{
					for (int k = 0; k < dimz; k++)
					{
						ind = i * dimy * dimz + j * dimz + k;
						nodes[ind].type = NODE_OUT;
					}
				}
				else
				{
					// set up & bottom bounds
					nodes[i * dimy * dimz + j * dimz + 0].type = NODE_OUT;
					for (int k = active_dimz-1; k < dimz; k++)
						nodes[i * dimy * dimz + j * dimz + k].type = NODE_OUT;

					nodes[i * dimy * dimz + j * dimz + 1].SetBound(BC_NOSLIP, BC_FREE, Vec3D(0.0, 0.0, 0.0), (FTYPE)baseT);
					nodes[i * dimy * dimz + j * dimz + active_dimz-2].SetBound(BC_NOSLIP, BC_FREE, Vec3D(0.0, 0.0, 0.0), (FTYPE)baseT);
					
					for (int k = 2; k < active_dimz-2; k++)
					{
						ind = i * dimy * dimz + j * dimz + k;
						switch (grid2D->GetType(i, j))
						{
						case NODE_BOUND:
							nodes[ind].SetBound(BC_NOSLIP, BC_FREE, Vec3D((FTYPE)grid2D->GetData(i, j).vel.x, (FTYPE)grid2D->GetData(i, j).vel.y, 0.0), (FTYPE)grid2D->GetData(i, j).T);
							break;
						case NODE_VALVE:
							if( grid2D->GetData(i, j).vel.x == 0 && grid2D->GetData(i, j).vel.y == 0 )
								nodes[ind].SetBound(BC_FREE, BC_FREE, Vec3D((FTYPE)grid2D->GetData(i, j).vel.x, (FTYPE)grid2D->GetData(i, j).vel.y, 0.0), (FTYPE)grid2D->GetData(i, j).T);
							else
								nodes[ind].SetBound(BC_NOSLIP, BC_NOSLIP, Vec3D((FTYPE)grid2D->GetData(i, j).vel.x, (FTYPE)grid2D->GetData(i, j).vel.y, 0.0), (FTYPE)grid2D->GetData(i, j).T);
							nodes[ind].type = NODE_VALVE;
							break;
						case NODE_IN:
							nodes[ind].type = NODE_IN;
							nodes[ind].T = (FTYPE)baseT;
							break;
						}
					}
				}
			}
	}

	double Grid3D::GetBarycentric(Vec2D v1, Vec2D v2, Vec2D v0, Vec2D p)
	{
		Vec2D en(v1.y - v2.y, v2.x - v1.x);
		double ec = (v1.x*v2.y - v1.y*v2.x);

		double fp = en.dot(p) + ec;
		double f0 = en.dot(v0) + ec;

		return fp/f0;
	}

	Vec2D GetIntersectHorizon(Vec2D p1, Vec2D p2, Vec2D p)
	{
		Vec2D res;
		res.y = p.y;
		if (fabs(p1.y - p2.y) < COMP_EPS) 
			res.x = p.x;
		else
			res.x = p1.x + (p2.x - p1.x) * (res.y - p1.y) / (p2.y - p1.y);
		return res;
	}

	// project the point back on the original 3D polygon and fill nodes
	void Grid3D::ProjectPointOnPolygon(DirType dir, int i, int j, Vec2D testp, Vec3D n, FTYPE d)
	{
		int k;
		switch( dir )
		{
		case X: 
			k = (int)((-d-testp.dot(Vec2D(n.y, n.z)))/n.x); 
			if( k >= 0 && k < dimx ) SetType(k, i, j, NODE_BOUND);
			break;
		case Y: 
			k = (int)((-d-testp.dot(Vec2D(n.x, n.z)))/n.y); 
			if( k >= 0 && k < dimy ) SetType(i, k, j, NODE_BOUND);
			break;
		case Z: 
			k = (int)((-d-testp.dot(Vec2D(n.x, n.y)))/n.z); 
			if( k >= 0 && k < dimz ) SetType(i, j, k, NODE_BOUND);
			break;
		}
	}

	void Grid3D::RasterPolygon(Vec3D p1, Vec3D p2, Vec3D p3, Vec3D v1, Vec3D v2, Vec3D v3, NodeType color)
    {
		// if zero polygon then immediate exit
		if( p1.equal(p2) && p1.equal(p3) ) return; 

		// compute normal
		Vec3D n = (p2-p1).cross(p3-p1);
		n.normalize();

		// compute plane distance
		FTYPE d = -p1.dot(n);

		// get max coordinate of the normal
		DirType dir;
		FTYPE maxv = max( fabs(n.x), max( fabs(n.y), fabs(n.z) ) );
		if( fabs(maxv-fabs(n.x)) < COMP_EPS ) dir = X;
		if( fabs(maxv-fabs(n.y)) < COMP_EPS ) dir = Y;
		if( fabs(maxv-fabs(n.z)) < COMP_EPS ) dir = Z;

		// project a polygon on the axis plane
		Vec2D pp1, pp2, pp3;
		int pdimx, pdimy;
		switch( dir )
		{
		case X: pp1 = Vec2D(p1.y, p1.z); pp2 = Vec2D(p2.y, p2.z); pp3 = Vec2D(p3.y, p3.z); pdimx = dimy; pdimy = dimz; break;
		case Y: pp1 = Vec2D(p1.x, p1.z); pp2 = Vec2D(p2.x, p2.z); pp3 = Vec2D(p3.x, p3.z); pdimx = dimx; pdimy = dimz; break;
		case Z: pp1 = Vec2D(p1.x, p1.y); pp2 = Vec2D(p2.x, p2.y); pp3 = Vec2D(p3.x, p3.y); pdimx = dimx; pdimy = dimy; break;
		}

		// sort points
		Vec2D mid;
		if( pp3.y < pp2.y ) { mid = pp3; pp3 = pp2; pp2 = mid; }
		if( pp1.y > pp2.y ) { mid = pp1; pp1 = pp2; pp2 = mid; }
		if( pp3.y < pp2.y ) { mid = pp3; pp3 = pp2; pp2 = mid; }
		
		// rasterize 2D triangle
		
		// compute mid point
		mid = GetIntersectHorizon(pp1, pp3, pp2);

		// compute slopes and steps
		Vec2D dir1(mid.x - pp1.x, mid.y - pp1.y);
		Vec2D dir2(pp3.x - mid.x, pp3.y - mid.y);
		int steps1 = (int)max(abs(dir1.x), abs(dir1.y)) + 1;
		int steps2 = (int)max(abs(dir2.x), abs(dir2.y)) + 1;
        Vec2D dp1((dir1.x) / steps1, (dir1.y) / steps1);		
		Vec2D dp2((dir2.x) / steps2, (dir2.y) / steps2);		
		
		Vec2D p = pp1;
		int di = (mid.x < pp2.x) ? 1 : -1;

		// go through the segment (pp1 - mid)
        for ( ; p.y < mid.y; ) 
        {
            int i = (int)p.x;
            int j = (int)p.y;

			int last_i = (int)GetIntersectHorizon(pp1, pp2, p).x;

			for (int i = (int)p.x; i != last_i + di; i += di)
				ProjectPointOnPolygon(dir, i, j, Vec2D((FTYPE)i, p.y), n, d);

			p += dp1;			
        }

		// go through the segment (mid - pp3)
        for ( ; p.y < pp3.y; ) 
        {
            int i = (int)p.x;
            int j = (int)p.y;

			int last_i = (int)GetIntersectHorizon(pp2, pp3, p).x;

			for (int i = (int)p.x; i != last_i + di; i += di)
				ProjectPointOnPolygon(dir, i, j, Vec2D((FTYPE)i, p.y), n, d);

			p += dp2;			
        }
	}

	void Grid3D::RasterLine(Vec3D p1, Vec3D p2, Vec3D v1, Vec3D v2, NodeType color)
    {
		Vec3D dir = p2 - p1;
		int steps = (int)max(abs(dir.x), max(abs(dir.y), abs(dir.z))) + 1;
        Vec3D dp = dir / (FTYPE)steps;
		
		Vec3D p = p1;
		
		// go through the whole segment
        for (int i = 0; i <= steps; i++) 
        {
            int x = (int)p.x;
            int y = (int)p.y;
			int z = (int)p.z;
			
			SetType(x, y, z, color);

			p += dp;
        }
	}

	void Grid3D::FloodFill(int start[3], NodeType color, int n, const int *neighborPos)
    {
		int *queue = new int[dimx * dimy * dimz * 3];
		int cur = -1;

		// we know that this cell is of our color type
		int last = 0;
		queue[0] = start[0];
		queue[1] = start[1];
		queue[2] = start[2];
		SetType(start[0], start[1], start[2], color);

		// run wave
		while (cur < last)
		{
			// get next index
			cur++;
			int i = queue[cur * 3 + 0];
			int j = queue[cur * 3 + 1];
			int k = queue[cur * 3 + 2];

			// add neighbors
			for (int t = 0; t < n; t++)
			{
				int next_i = i + neighborPos[t * 3 + 0];
				int next_j = j + neighborPos[t * 3 + 1];
				int next_k = k + neighborPos[t * 3 + 2];

				if ((next_i >= 0) && (next_i < dimx) && (next_j >= 0) && (next_j < dimy) && (next_k >= 0) && (next_k < dimz))
					if (GetType(next_i, next_j, next_k) == NODE_IN) 
					{
						last++;
						queue[last * 3 + 0] = next_i;
						queue[last * 3 + 1] = next_j;
						queue[last * 3 + 2] = next_k;
						SetType(next_i, next_j, next_k, color);
					}
			}
		}

		delete [] queue;
	}

	void Grid3D::Build(FrameInfo3D &frame)
	{
		// mark all cells as inner 
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				for (int k = 0; k < dimz; k++)
					SetType(i, j, k, NODE_IN);
     
		// rasterize polygons (boundary)
		for( int s = 0; s < frame.NumShapes; s++ )
			for( int i = 0; i < frame.Shapes[s].NumIndices; i++ )
				if( !frame.Shapes[s].Active )
				{
					int i1 = frame.Shapes[s].Indices[i*3+0];
					int i2 = frame.Shapes[s].Indices[i*3+1];
					int i3 = frame.Shapes[s].Indices[i*3+2];
					RasterPolygon(frame.Shapes[s].Vertices[i1], frame.Shapes[s].Vertices[i2], frame.Shapes[s].Vertices[i3],
								  frame.Shapes[s].Velocities[i1], frame.Shapes[s].Velocities[i2], frame.Shapes[s].Velocities[i3],
								  NODE_BOUND);

					// additionally rasterize all edges to cover holes
					RasterLine(frame.Shapes[s].Vertices[i1], frame.Shapes[s].Vertices[i2], frame.Shapes[s].Velocities[i1], frame.Shapes[s].Velocities[i2], NODE_BOUND);
					RasterLine(frame.Shapes[s].Vertices[i1], frame.Shapes[s].Vertices[i3], frame.Shapes[s].Velocities[i1], frame.Shapes[s].Velocities[i3], NODE_BOUND);
					RasterLine(frame.Shapes[s].Vertices[i3], frame.Shapes[s].Vertices[i2], frame.Shapes[s].Velocities[i3], frame.Shapes[s].Velocities[i2], NODE_BOUND);
				}

		// detect all outside nodes by running wave algorithm
		const int neighborPos[18] = { -1, 0, 0,  1, 0, 0,  0, -1, 0,  0, 1, 0,  0, 0, -1,  0, 0, 1 };
		int start[3] = { 0, 0, 0 };
		FloodFill(start, NODE_OUT, 6, neighborPos); 
		
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				for (int k = 0; k < dimz; k++)
				{
					NodeType t = GetType(i, j, k);
					switch (t)
					{	
						case NODE_IN: 
						case NODE_OUT: 
							SetData(i, j, k, BC_NOSLIP, BC_NOSLIP, Vec3D(0, 0, 0), (FTYPE)baseT); 
							break; 
					}
				}
	}

	void Grid3D::ComputeSubframeInfo(int frame, FTYPE substep, FrameInfo3D &res)
	{
		int framep1 = (frame + 1) % num_frames;

		res.Duration = 0;

		FTYPE isubstep = 1 - substep;

		res.Init(frames[frame].NumShapes);
		for (int i=0; i<res.NumShapes; i++)
		{
			res.Shapes[i].Init(frames[frame].Shapes[i]);
			for (int k=0; k<res.Shapes[i].NumVertices; k++)
			{
				res.Shapes[i].Vertices[k] = frames[frame].Shapes[i].Vertices[k] * isubstep + frames[framep1].Shapes[i].Vertices[k] * substep;
				res.Shapes[i].Velocities[k] = frames[frame].Shapes[i].Velocities[k] * isubstep + frames[framep1].Shapes[i].Velocities[k] * substep;
				res.Shapes[i].Active = frames[frame].Shapes[i].Active;
			}

			for (int k=0; k<res.Shapes[i].NumIndices; k++)
				for (int j=0; j<3; j++)
					res.Shapes[i].Indices[k*3+j] = frames[frame].Shapes[i].Indices[k*3+j];
		}

		if (frames[frame].Field.Correlate(frames[framep1].Field))
		{
			res.Field.Init(frames[frame].Field);
			int nx = res.Field.Nx;
			int ny = res.Field.Ny;
			int nz = res.Field.Nz;

			for (int k=0; k<nz; k++)
				for (int j=0; j<ny; j++)
					for (int i=0; i<nx; i++)
					{
						int t = k * nx * ny + j * nx + i;

						Vec3D v1 = frames[frame].Field.Data[t];
						Vec3D v2 = frames[framep1].Field.Data[t];

						res.Field.Data[t] = Vec3D(0, 0, 0);
						if ( ( fabs(v1.length() ) > COMP_EPS ) && ( fabs(v2.length() ) > COMP_EPS ) )
							res.Field.Data[t] = v1 * isubstep + v2 * substep;						
					};
		}
	}

	void Grid3D::Prepare3D_Shape(double time)
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

		FrameInfo3D info;
		ComputeSubframeInfo(frame % num_frames, (FTYPE)substep, info);
		Build(info);
	}

	void Grid3D::Prepare3D_NetCDF(double time)
	{
		// mark all cells as inner 
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
				for (int k = 0; k < dimz; k++)
				{
					SetType(i, j, k, NODE_OUT);
					SetData(i, j, k, BC_NOSLIP, BC_NOSLIP, Vec3D(0, 0, 0), (FTYPE)baseT); 
				}
     
		// mark sea cells (depth < 0)
		for (int i = 0; i < dimx; i++)
			for (int j = 0; j < dimy; j++)
			{
				// compute corresponding input location
				int di = i * depthInfo->dimx / dimx;
				int dj = j * depthInfo->dimy / dimy;

				// compute sea depth
				float z = depthInfo->depth[dj + di * depthInfo->dimy];

				if (z < 0.0f)
				{
					int bound_k = (int)(dimz * z / bbox.pMin.z);
					for (int k = 1; k < bound_k; k++)
					{
						int ind = i * dimy * dimz + j * dimz + k;
						nodes[ind].type = NODE_IN;
					}
				}
			}

		// compute boundaries
		for (int i = 1; i < dimx-1; i++)
			for (int j = 1; j < dimy-1; j++)
				for (int k = 1; k < dimz-1; k++)
					if (GetType(i, j, k) == NODE_IN)
					{
						bool bound = (GetType(i-1, j, k) == NODE_OUT) || (GetType(i+1, j, k) == NODE_OUT) ||
									 (GetType(i, j-1, k) == NODE_OUT) || (GetType(i, j+1, k) == NODE_OUT) ||
									 (GetType(i, j, k-1) == NODE_OUT) || (GetType(i, j, k+1) == NODE_OUT);
						if( bound ) 
						{
							int ind = i * dimy * dimz + j * dimz + k;
							nodes[ind].SetBound(BC_NOSLIP, BC_NOSLIP, Vec3D(0.0f, 0.0f, 0.0f), (float)baseT);
						}
					}

		int num = 0;
		int *indices = new int[(dimx * dimy + dimy * dimz + dimz * dimx) * 2];

		for (int i = 1; i < dimx-1; i++)
			for (int j = 1; j < dimy-1; j++)
				for (int k = 1; k < dimz-1; k++)
					if (GetType(i, j, k) == NODE_OUT)
					{
						bool bound = (GetType(i-1, j, k) == NODE_BOUND) || (GetType(i+1, j, k) == NODE_BOUND) ||
									 (GetType(i, j-1, k) == NODE_BOUND) || (GetType(i, j+1, k) == NODE_BOUND) ||
									 (GetType(i, j, k-1) == NODE_BOUND) || (GetType(i, j, k+1) == NODE_BOUND);
						if( bound ) 
						{
							int ind = i * dimy * dimz + j * dimz + k;
							indices[num++] = ind;
						}
					}

		for (int i = 0; i < num; i++)
			nodes[indices[i]].SetBound(BC_NOSLIP, BC_NOSLIP, Vec3D(0.0f, 0.0f, 0.0f), (float)baseT);
		delete [] indices;

		// set input/output streams on quad boundaries (special for sea test)
		// currently upper stream is in, lower stream is out
		int start, end;
		for (int i = 0; i < dimx; i++) {
			start = -1;
			for (int k = 0; k < dimz; k++) 
				if( GetType(i, dimy-1, k) == NODE_IN ) {
					if( start < 0 ) start = k;
					end = k;
				}
			for (int k = 0; k < dimz; k++)
				if (GetType(i, dimy-1, k) == NODE_IN) 
				{
					SetType(i, dimy-1, k, NODE_VALVE);
					if( k < (start+end)/2 )
						SetData(i, dimy-1, k, BC_NOSLIP, BC_NOSLIP, bcInVel, (float)bcInT);
					else
						SetData(i, dimy-1, k, BC_NOSLIP, BC_NOSLIP, Vec3D()-bcInVel, 2.0f-(float)bcInT);
				}
		}

		for (int j = 0; j < dimy; j++) {
			start = -1;
			for (int k = 0; k < dimz; k++) 
				if( GetType(dimx-1, j, k) == NODE_IN ) {
					if( start < 0 ) start = k;
					end = k;
				}
			for (int k = 0; k < dimz; k++)
				if (GetType(dimx-1, j, k) == NODE_IN) 
				{
					SetType(dimx-1, j, k, NODE_VALVE);
					if( k < (start+end)/2 )
						SetData(dimx-1, j, k, BC_NOSLIP, BC_NOSLIP, bcInVel, (float)bcInT);
					else
						SetData(dimx-1, j, k, BC_NOSLIP, BC_NOSLIP, Vec3D()-bcInVel, 2.0f-(float)bcInT);
				}
		}
	}

	void Grid3D::TestPrint(char *filename)
	{
		PARAplan *pplan = PARAplan::Instance();
		if (pplan->rank() > 0)
			return;
		FILE *file = NULL;
		fopen_s(&file, filename, "w");
		fprintf(file, "grid (z-slices):\n");
		fprintf(file, "%i %i %i\n", dimx, dimy, dimz);
		for (int k = 0; k < dimz; k++)
		{
			fprintf(file, "%i\n", k);
			for (int i = 0; i < dimx; i++)
			{
				for (int j = 0; j < dimy; j++)
				{					
					if( k == 55 && j == 2 && i == 0 ) { 
						fprintf(file, "!"); 
						continue;
					}
					NodeType t = GetType(i, j, k);
					switch (t)
					{
					case NODE_IN: fprintf(file, " "); break;
					case NODE_OUT: fprintf(file, "."); break;
					case NODE_BOUND: fprintf(file, "#"); break;		
					case NODE_VALVE: fprintf(file, "+"); break;		
					}
				}
				fprintf(file, "\n");
			}
		}
		fclose(file);
	}

	void Grid3D::OutputImage(char *filename_base)
	{
		PARAplan *pplan = PARAplan::Instance();
		if (pplan->rank() > 0)
			return;
		BitmapFileHeader bfh;
		BitmapInfoHeader bih;
		
		memset(&bfh, sizeof(bfh), 0);
		bfh.bfType = 0x4D42;
		bfh.bfOffBits = sizeof(bfh) + sizeof(bih);
		bfh.bfSize = bfh.bfOffBits + sizeof(char) * 3 * dimx * dimy + dimx * (dimy % 4);
		
		memset(&bih, sizeof(bih), 0);
		bih.biSize = int(sizeof(bih));
		bih.biBitCount = 24;
		bih.biCompression = 0;
		bih.biHeight = dimx;
		bih.biWidth = dimy;
		bih.biPlanes = 1;

		//bih.biSizeImage = bfh.bfSize - bfh.bfOffBits;
		bih.biClrUsed = 8;

		for (int k = 0; k < dimz; k++)
		{
			FILE *file = NULL;
			string filename = filename_base;

			filename.append("/");
			filename.append(stringify(k));
			filename.append(".bmp");
			if( fopen_s(&file, filename.c_str(), "w") )
			{
				_mkdir(filename_base);
				fopen_s(&file, filename.c_str(), "w");
			}
			
			fwrite(&bfh, sizeof(bfh), 1, file);
			fwrite(&bih, sizeof(bih), 1, file);

			for (int i = dimx-1; i >= 0; i--)
			{
				char color[3];
				for (int j = 0; j < dimy; j++)
				{
					NodeType t = GetType(i, j, k);
					switch (t)
					{
					case NODE_IN: color[0] = (char)245; color[1] = (char)73; color[2] = (char)69; break;		// blue
					case NODE_OUT: color[0] = (char)0; color[1] = (char)0; color[2] = (char)0; break;			// black
					case NODE_BOUND: color[0] = (char)255; color[1] = (char)255; color[2] = (char)255; break;	// white
					case NODE_VALVE: color[0] = (char)241; color[1] = (char)41; color[2] = (char)212; break;	// purple
					}
					fwrite(color, sizeof(char) * 3, 1, file);
				}
				fwrite(color, sizeof(char) * 3, dimy % 4, file);
			}

			fclose(file);
		}
	}
}
