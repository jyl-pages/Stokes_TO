#include <iostream>
#include "MeshAdvFunc.h"
#include "MeshFunc.h"
#include "Mesh.h"
#include "Constants.h"
#include "AuxFunc.h"

#ifdef USE_TRI2D
#include "Triangulation2D.h"
#endif
#ifdef USE_TETGEN
#include "TetGenFunc.h"
#endif

namespace MeshFunc{
	void Initialize_Circle_Mesh(const Vector2& c,const real r,TriangleMesh<2>* mesh,const int sub_n)
	{
		#ifdef USE_TRI2D
		Array<Array<Vector2> > loops(1);
		real theta=(real)two_pi/(real)sub_n;
		real max_seg_length=two_pi*r/(real)sub_n;
		for(int i=0;i<sub_n;i++){
			real t=theta*(real)i;
			Vector2 pos(r*cos(t),r*sin(t));loops[0].push_back(c+pos);}
		Triangulation::Triangulation2D(loops,*mesh,max_seg_length);
		#else
		std::cerr << "Error: [MeshFunc] Initialize_Circle_Mesh requires USE_TRI2D" << std::endl;
		#endif
	}

	void Initialize_Circle_Mesh(const Vector3& c,const real r,const Vector3& normal,TriangleMesh<3>* mesh,const int sub_n)
	{
		TriangleMesh<2> mesh_2d;
		Initialize_Circle_Mesh(Vector2::Zero(),r,&mesh_2d,sub_n);
		int vn=(int)mesh_2d.Vertices().size();
		mesh->Vertices().resize(vn);
		mesh->Elements()=mesh_2d.Elements();
		Vector3 v1=Vector3::Unit(2);
		Vector3 axis;real theta;AuxFunc::Angle_And_Axis_Between(v1,normal,theta,axis);
		AngleAxis rot(theta,axis);

		for(int i=0;i<vn;i++){
			const Vector2& v=mesh_2d.Vertices()[i];
			mesh->Vertices()[i]=rot*Vector3(v[0],v[1],(real)0)+c;}
	}

	void Initialize_Ring_Mesh(const Vector2& c, const real R, const real r, TriangleMesh<2>* mesh, const int sub_n)
	{
		#ifdef USE_TRI2D
		Array<Array<Vector2> > loops(2);
		real theta=(real)two_pi/(real)sub_n;
		real max_seg_length=two_pi*r/(real)sub_n;
		for(int i=0;i<sub_n;i++){
			real t=theta*(real)i;
			Vector2 pos_r(r*cos(t),r*sin(t));
			Vector2 pos_R(R*cos(t),R*sin(t));
			loops[0].push_back(c+pos_R);
			loops[1].push_back(c+pos_r);}
		Triangulation::Triangulation2D(loops,*mesh,max_seg_length);
		#else
		std::cerr << "Error: [MeshFunc] Initialize_Circle_Mesh requires USE_TRI2D" << std::endl;
		#endif
	}

	void Initialize_Polygon_Mesh(Array<Vector2>& vtx,TriangleMesh<2>* mesh,const real dx)
	{
		#ifdef USE_TRI2D
		Array<Array<Vector2> > loops(1);int n=(int)vtx.size();
		for(int i=0;i<n;i++){
			Vector2 dir=vtx[(i+1)%n]-vtx[i];
			real length=dir.norm();dir/=length;
			int m=(int)(length/dx);real step=length/(real)m;
			for(int j=0;j<m;j++){ loops[0].push_back(vtx[i]+dir*(real)j*step);}}
		Triangulation::Triangulation2D(loops,*mesh,(real)1.2*dx);
		#else
		std::cerr << "Error: [MeshFunc] Initialize_Polygon_Mesh requires USE_TRI2D" << std::endl;
		#endif
	}

	void Initialize_Polygon_Mesh(const Array<Array<Vector2> >& vtx,TriangleMesh<2>* mesh,const real dx)
	{
		#ifdef USE_TRI2D
		Array<Array<Vector2> > loops(vtx.size());
		for(int p=0;p<loops.size();p++){
			int n=(int)vtx[p].size();
			for(int i=0;i<n;i++){
				Vector2 dir=vtx[p][(i+1)%n]-vtx[p][i];
				real length=dir.norm();dir/=length;
				int m=(int)(length/dx);real step=length/(real)m;
				for(int j=0;j<m;j++){loops[p].push_back(vtx[p][i]+dir*(real)j*step);}}}

		Triangulation::Triangulation2D(loops,*mesh,(real)1.2*dx);
		#else
		std::cerr << "Error: [MeshFunc] Initialize_Polygon_Mesh requires USE_TRI2D" << std::endl;
		#endif
	}

	void Initialize_Polygon_Mesh_From_File(const char* file_name, TriangleMesh<2>* mesh) {
		Array<Array<Vector2>> vtx;
		real dx = 0.05;
		if (!Read_Vtx_From_File(file_name, vtx, 2, dx)) { return; }
		std::cout << "Finished reading polygon file with " << vtx.size() << " loop(s)" << std::endl;
		Initialize_Polygon_Mesh(vtx, mesh, dx);
	}

	void Initialize_Polygon_Mesh_From_File(const char* file_name, TetrahedronMesh<3>* mesh) {
		Array<Array<Vector3>> vtx;
		real dx = 0.05;
		if (!Read_Vtx_From_File(file_name, vtx, 3, dx)) { return; }
		std::cout << "Finished reading polygon file with " << vtx.size() << "loop(s)" << std::endl;
		Initialize_Polygon3d_Mesh(vtx, *mesh, dx);
	}

	template<typename T> bool Read_Vtx_From_File(const char* file_name, Array<Array<T>>& vtx, int d, real& dx, int mode) {
		std::ifstream in(file_name); char line[128] = { 0 };
		int recentLoop = -1; int recentVtx = -1;
		while (in.getline(line, sizeof(line))) {
			std::stringstream content(line);
			if (recentLoop == -1) {
				std::string head = ""; int loopNum;
				if (mode == 0) {
					if (content >> head >> loopNum >> dx && !head.compare(char(d + '0') + std::string("dContour"))) {
						std::cout << "Reading 2D contour file with " << loopNum << " loop(s)!" << std::endl;
						vtx = Array<Array<T>>(loopNum); recentLoop++;
					}
					else { std::cerr << "Invalid file head: " << head << std::endl; return false; }
				}
				else if (mode == 1) {
					int vtxNum;
					if (content >> head >> vtxNum && !head.compare(char(d + '0') + std::string("dPoints"))) {
						std::cout << "Reading 2D point sets file with " << vtxNum << " points!" << std::endl;
						loopNum = 1;
						vtx = Array<Array<T>>(loopNum); recentLoop++; vtx[recentLoop] = Array<T>(vtxNum); recentVtx++;
					}
					else { std::cerr << "Invalid file head: " << head << std::endl; return false; }
				}

			}
			else if (recentVtx == -1) {
				int vtxNum;
				if (content >> vtxNum) { vtx[recentLoop] = Array<T>(vtxNum); recentVtx++; }
				else { std::cerr << "Invalid contour!" << std::endl; return false; }
			}
			else {
				T x = T();
				if (content >> x[0] >> x[1] && (d == 2 || content >> x[2])) {
					vtx[recentLoop][recentVtx++] = x;
					if (recentVtx == vtx[recentLoop].size()) { recentLoop++; recentVtx = -1; }
				}
				else { std::cerr << "Invalid vertex!" << std::endl; return false; }
			}
		}
		if (recentLoop == -1) { std::cerr << "Reading failed!" << std::endl; return false; }
		std::cout << "Finished reading vertices" << std::endl; in.close();
		return true;
	}

	void Initialize_Rectangle_Mesh(const Vector2& domain_min,const Vector2& domain_size,TriangleMesh<2>* mesh,const real dx)
	{
		#ifdef USE_TRI2D
		real w=domain_size[0];real h=domain_size[1];
		Array<Vector2> vtx={domain_min,domain_min+Vector2::Unit(0)*w,domain_min+Vector2::Unit(0)*w+Vector2::Unit(1)*h,domain_min+Vector2::Unit(1)*h};
		Initialize_Polygon_Mesh(vtx,mesh,dx);
		#else
		std::cerr << "Error: [MeshFunc] Initialize_Rectangle_Mesh requires USE_TRI2D" << std::endl;
		#endif
	}

	void Initialize_Sphere_Mesh(const Vector3& c, const real r, TetrahedronMesh<3>& volume_mesh, const real dx)
	{
		#ifdef USE_TETGEN
		TriangleMesh<3> surface_mesh;MeshFunc::Initialize_Sphere_Mesh(r,&surface_mesh,1);
		TetGen::Generate(surface_mesh,volume_mesh,dx);
		for(int i=0;i<volume_mesh.Vertices().size();i++)volume_mesh.Vertices()[i]+=c;
		#else
		std::cerr<<"Error: [MeshFunc] Initialize_Sphere_Mesh requires USE_TETGEN"<<std::endl;
		#endif
	}

	void Initialize_Cuboid_Mesh(const Vector3& domain_min, const Vector3& domain_size, TetrahedronMesh<3>& volume_mesh, const real dx)
	{
		#ifdef USE_TETGEN
		real w = domain_size[0]; real h = domain_size[1]; real l = domain_size[2];
		Array<Vector3> vtx = { domain_min,domain_min + Vector3::Unit(0) * w,domain_min + Vector3::Unit(0) * w + Vector3::Unit(1) * h,domain_min + Vector3::Unit(1) * h,
		domain_min + Vector3::Unit(2) * l,domain_min + Vector3::Unit(2) * l + Vector3::Unit(0) * w,domain_min + Vector3::Unit(2) * l + Vector3::Unit(0) * w + Vector3::Unit(1) * h,domain_min + Vector3::Unit(2) * l + Vector3::Unit(1) * h };
		int ele[12][3] = { {1,2,4},{2,3,4},{3,4,7},{4,7,8},{2,3,7},{2,6,7},{1,4,5},{4,5,8},{1,2,5},{2,5,6},{5,6,7},{5,7,8} };
		TriangleMesh<3> surface_mesh; surface_mesh.Vertices().resize(8); surface_mesh.Elements().resize(12);
		for (int i = 0; i < 8; i++)surface_mesh.Vertices()[i] = vtx[i];
		for (int i = 0; i < 12; i++)surface_mesh.Elements()[i] = Vector3i(ele[i][0]-1, ele[i][1]-1, ele[i][2]-1);
		TetGen::Generate(surface_mesh, volume_mesh, dx);
		#else
		std::cerr << "Error: [MeshFunc] Initialize_Cuboid_Mesh requires USE_TETGEN" << std::endl;
		#endif
	}

	void Initialize_Cylinder_Mesh(const Vector3& c, const real r, const real h, TetrahedronMesh<3>& volume_mesh, const real dx)
	{
		//c:center of cylinder, r： radius of cylinder, h: height of cylinder
		#ifdef USE_TETGEN
		Array<Vector3> vtx;
		vtx.push_back(c + Vector3(0, 0, h / 2));
		for (int i = 0; i < 12; i++)
		{
			vtx.push_back(c + Vector3(r * cos(i * pi / 6), r * sin(i * pi / 6), h / 2));
		}
		vtx.push_back(c + Vector3(0, 0, -h / 2));
		for (int i = 0; i < 12; i++)
		{
			vtx.push_back(c + Vector3(r * cos(i * pi / 6), r * sin(i * pi / 6), -h / 2));
		}
		TriangleMesh<3> surface_mesh; surface_mesh.Vertices().resize(26); surface_mesh.Elements().resize(48);
		int ele[48][3] = { {1,2,3},{1,3,4},{1,4,5},{1,5,6},{1,6,7},{1,7,8},{1,8,9},{1,9,10},{1,10,11},{1,11,12},{1,12,13},{1,13,2},
		{14,15,16},{14,16,17},{14,17,18},{14,18,19},{14,19,20},{14,20,21},{14,21,22},{14,22,23},{14,23,24},{14,24,25},{14,25,26},{14,26,15},
		{2,3,15},{3,15,16},{3,4,16},{4,16,17},{4,5,17},{5,17,18},{5,6,18},{6,18,19}, {6,7,19},{7,19,20},{7,8,20},{8,20,21},{8,9,21},{9,21,22},
		{9,10,22},{10,22,23},{10,11,23},{11,23,24},{11,12,24},{12,24,25},{12,13,25},{13,25,26},{13,2,26},{2,15,26} };
		for (int i = 0; i < 26; i++) surface_mesh.Vertices()[i] = vtx[i];
		for (int i = 0; i < 48; i++)surface_mesh.Elements()[i] = Vector3i(ele[i][0] - 1, ele[i][1] - 1, ele[i][2] - 1);
		TetGen::Generate(surface_mesh, volume_mesh, dx);
		#else
		std::cerr << "Error: [MeshFunc] Initialize_Cylinder_Mesh requires USE_TETGEN" << std::endl;
		#endif
	}

	void Initialize_Polygon3d_Mesh(const Array<Array<Vector3>>& vtx, TetrahedronMesh<3>& volume_mesh, const real dx)
	{
		#ifdef USE_TETGEN
		TriangleMesh<3> surface_mesh;

		TetGen::Generate(surface_mesh, volume_mesh, dx);
		#else
		std::cerr << "Error: [MeshFunc] Initialize_Polygon3d_Mesh requires USE_TETGEN" << std::endl;
		#endif
	}
};