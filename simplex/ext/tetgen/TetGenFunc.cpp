//#####################################################################
// TetGen Func
// Jingyuan Zhu
//#####################################################################
#include <iostream>
#include <memory>
#include "tetgen.h"
#include "TetGenFunc.h"
#include <cstdlib>
#include <fstream>

namespace TetGen{

	void Generate(const TriangleMesh<3>& surface_mesh, TetrahedronMesh<3>& volume_mesh, const real dx)
	{
		//Generate Input File
		std::fstream fp;
		fp.open("data.poly", std::ios::out);
		fp << (int)surface_mesh.Vertices().size() << "  " << (int)3 << std::endl;
		for (int i = 0; i < surface_mesh.Vertices().size(); i++)
		{
			fp << i << " " << surface_mesh.Vertices()[i].x() << "  " << surface_mesh.Vertices()[i].y() << "  " << surface_mesh.Vertices()[i].z() << std::endl;
		}
		fp << std::endl;
		fp << (int)surface_mesh.Elements().size() << "  " << (int)0 << std::endl;
		for (int i = 0; i < surface_mesh.Elements().size(); i++)
		{
			fp << (int)1 << "  " << (int)0 << "  " << (int)0 << std::endl;
			fp << (int)3 << "  " << surface_mesh.Elements()[i][0] << "  " << surface_mesh.Elements()[i][1] << "  " << surface_mesh.Elements()[i][2] << std::endl;
		}
		//hole and region
		fp << (int)0 << std::endl;
		fp << (int)0 << std::endl;
		fp.close();
		//Write mtr File
		if (dx != 0)
		{
			fp.open("data.mtr", std::ios::out);
			fp << surface_mesh.Vertices().size() << "  " << (int)1 << std::endl;
			for (int i = 0; i < surface_mesh.Vertices().size(); i++)
			{
				fp << dx << std::endl;
			}
			fp.close();
		}
		//Tetgen generate TetMesh
		tetgenbehavior b;
		b.plc = 1;	//-p
		tetgenio in, addin, bgmin;
		strcpy(b.infilename, "data");
		in.load_plc(b.infilename, 1);
		if (dx != 0)
		{
			b.quality = 1; //Refine Mesh -q
			b.metric = 1;//-m
			//strcpy(b.bgmeshfilename, "data.b");
			//bgmin.load_tetmesh(b.bgmeshfilename, 1);
		}
		tetrahedralize(&b, &in, NULL, &addin, &bgmin);
		//Read TetMesh from Tetgen OutputFile
		std::ifstream nodefile, elefile;
		nodefile.open(".node", std::ios::in);
		if (!nodefile.is_open())
			std::cout << "Open Node File Failure" << std::endl;
		elefile.open(".ele", std::ios::in);
		if (!elefile.is_open())
			std::cout << "Open Ele File Failure" << std::endl;
		Array<Vector3>& _vertices = volume_mesh.Vertices();
		Array<Vector4i>& _elements = volume_mesh.Elements();
		int node_size, node_dim, node_index, node_hole, node_prop;
		nodefile >> node_size >> node_dim >> node_hole >> node_prop;
		_vertices.resize(node_size);
		for (int i = 0; i < node_size; i++)
			nodefile >> node_index >> _vertices[i][0] >> _vertices[i][1] >> _vertices[i][2];
		/*for (int i = 0; i < node_size; i++)
			std::cout << _vertices[i][0] << "  " << _vertices[i][1] << "  " << _vertices[i][2] << std::endl;*/
		nodefile.close();
		int ele_size, ele_nodepertet, ele_prop, ele_index;
		elefile >> ele_size >> ele_nodepertet >> ele_prop;
		_elements.resize(ele_size);
		for (int i = 0; i < ele_size; i++)
			elefile >> ele_index >> _elements[i][0] >> _elements[i][1] >> _elements[i][2] >> _elements[i][3];
		elefile.close();
		/*for (int i = 0; i < ele_size; i++)
			std::cout << _elements[i][0] << "  " << _elements[i][1] << "  " << _elements[i][2] << "  " << _elements[i][3] << std::endl;*/
		remove("data.poly"); remove("data.mtr");
		remove(".ele"); remove(".node"); remove(".face");remove(".edge"); remove(".mtr"); remove(".p2t");
	}
};