//////////////////////////////////////////////////////////////////////////
// Marching cube
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "LevelSet.h"
#include "MarchingCubes.h"

MarchingCubes<2>::MarchingCubes(LevelSet<2>& _levelset,std::shared_ptr<SurfaceMesh<2> > _mesh/*=nullptr*/,const real _contour_value/*=(real)0*/)
:levelset(_levelset),grid(_levelset.grid),contour_value(_contour_value)
{
	if(_mesh==nullptr)mesh=std::make_shared<SurfaceMesh<2> >();else mesh=_mesh;
}

void MarchingCubes<2>::Marching()
{
	const VectorDi& cell_counts=grid.cell_counts-VectorDi::Ones();const VectorDi& node_counts=grid.cell_counts;
	Array<VectorD>& vertices=(*mesh->vertices);vertices.clear();
	Array<VectorDi>& edges=mesh->elements;edges.clear();
	Array<int> vtx_idx_left;vtx_idx_left.resize(cell_counts[1]);std::fill(vtx_idx_left.begin(),vtx_idx_left.end(),-1);
	Array<int> vtx_idx_right;vtx_idx_right.resize(cell_counts[1]);std::fill(vtx_idx_right.begin(),vtx_idx_right.end(),-1);
	Array<int> vtx_idx_updown;vtx_idx_updown.resize(node_counts[1]);std::fill(vtx_idx_updown.begin(),vtx_idx_updown.end(),-1);
	////calculate column 1 left
	for(int j=0;j<cell_counts[1];j++){
		real val_up=levelset.phi(VectorDi(0,j+1)),val_down=levelset.phi(VectorDi(0,j));
		if(Edge_Has_Vertex(val_up,val_down,contour_value)){
			real alpha=(val_up-contour_value)/(val_up-val_down);VectorD pos=((real)1-alpha)*grid.Center(VectorDi(0,j+1))+alpha*grid.Center(VectorDi(0,j));
			vertices.push_back(pos);vtx_idx_right[j]=(int)vertices.size()-1;}}
	////calculate other columns
	for(int i=0;i<cell_counts[0];i++){
		vtx_idx_left=vtx_idx_right;
		std::fill(vtx_idx_right.begin(),vtx_idx_right.end(),-1);
		std::fill(vtx_idx_updown.begin(),vtx_idx_updown.end(),-1);
		////calculate up,down vertices
		for(int j=0;j<node_counts[1];j++){
			real val_left=levelset.phi(VectorDi(i,j)),val_right=levelset.phi(VectorDi(i+1,j));
			if(Edge_Has_Vertex(val_left,val_right,contour_value)){
			real alpha=(val_left-contour_value)/(val_left-val_right);VectorD pos=(1-alpha)*grid.Center(VectorDi(i,j))+alpha*grid.Center(VectorDi(i+1,j));
			vertices.push_back(pos);vtx_idx_updown[j]=(int)vertices.size()-1;}}
		////calculate right vertices
		for(int j=0;j<cell_counts[1];j++){
			real val_up=levelset.phi(VectorDi(i+1,j+1)),val_down=levelset.phi(VectorDi(i+1,j));
			if(Edge_Has_Vertex(val_up,val_down,contour_value)){
				real alpha=(val_up-contour_value)/(val_up-val_down);VectorD pos=(1-alpha)*grid.Center(VectorDi(i+1,j+1))+alpha*grid.Center(VectorDi(i+1,j));
				vertices.push_back(pos);vtx_idx_right[j]=(int)vertices.size()-1;}}
		////calculate edges
		for(int j=0;j<cell_counts[1];j++){
			real value[4]; unsigned int square_type;
			value[0]=levelset.phi(VectorDi(i,j));value[1]=levelset.phi(VectorDi(i+1,j));
			value[2]=levelset.phi(VectorDi(i+1,j+1));value[3]=levelset.phi(VectorDi(i,j+1));
			square_type=Square_Type(value[0],value[1],value[2],value[3],contour_value);
			switch(square_type){
				case 0:
				case 15:/*do nothing*/
					break;
				case 1:edges.push_back(VectorDi(vtx_idx_updown[j],vtx_idx_left[j]));
					break;
				case 14:edges.push_back(VectorDi(vtx_idx_left[j],vtx_idx_updown[j]));
					break;
				case 2:edges.push_back(VectorDi(vtx_idx_right[j],vtx_idx_updown[j]));
					break;
				case 13:edges.push_back(VectorDi(vtx_idx_updown[j],vtx_idx_right[j]));
					break;
				case 3:edges.push_back(VectorDi(vtx_idx_right[j],vtx_idx_left[j]));
					break;
				case 12:edges.push_back(VectorDi(vtx_idx_left[j],vtx_idx_right[j]));
					break;
				case 4:edges.push_back(VectorDi(vtx_idx_updown[j+1],vtx_idx_right[j]));
					break;
				case 11:edges.push_back(VectorDi(vtx_idx_right[j],vtx_idx_updown[j+1]));
					break;
				case 5:
					{real avg_v=(real).25*(levelset.phi(VectorDi(i,j))+levelset.phi(VectorDi(i+1,j))+levelset.phi(VectorDi(i+1,j+1))+levelset.phi(VectorDi(i,j+1)));
					if(avg_v*levelset.phi(VectorDi(i,j))>contour_value){
						edges.push_back(VectorDi(vtx_idx_updown[j+1],vtx_idx_right[j]));
						edges.push_back(VectorDi(vtx_idx_updown[j],vtx_idx_left[j]));}
					else{
						edges.push_back(VectorDi(vtx_idx_updown[j+1],vtx_idx_left[j]));
						edges.push_back(VectorDi(vtx_idx_updown[j],vtx_idx_right[j]));}}
					break;
				case 10:
					{real avg_v=(real).25*(levelset.phi(VectorDi(i,j))+levelset.phi(VectorDi(i+1,j))+levelset.phi(VectorDi(i+1,j+1))+levelset.phi(VectorDi(i,j+1)));
					if(avg_v*levelset.phi(VectorDi(i,j))>contour_value){
						edges.push_back(VectorDi(vtx_idx_right[j],vtx_idx_updown[j+1]));
						edges.push_back(VectorDi(vtx_idx_left[j],vtx_idx_updown[j]));}
					else{
						edges.push_back(VectorDi(vtx_idx_left[j],vtx_idx_updown[j+1]));
						edges.push_back(VectorDi(vtx_idx_right[j],vtx_idx_updown[j]));}}
					break;
				case 6:edges.push_back(VectorDi(vtx_idx_updown[j+1],vtx_idx_updown[j]));
					break;
				case 9:edges.push_back(VectorDi(vtx_idx_updown[j],vtx_idx_updown[j+1]));
					break;
				case 7:edges.push_back(VectorDi(vtx_idx_updown[j+1],vtx_idx_left[j]));
					break;
				case 8:edges.push_back(VectorDi(vtx_idx_left[j],vtx_idx_updown[j+1]));
					break;
				default:std::cerr<<"Invalid marching square type!"<<std::endl;
					break;}}}
}

unsigned int MarchingCubes<2>::Square_Type(real v0,real v1,real v2,real v3,real iso_value/*=(real)0*/)
{unsigned int type=0;if(v0<iso_value) type|=1;if(v1<iso_value) type|=2;if(v2<iso_value) type|=4;if(v3<iso_value) type|=8;return type;}

bool MarchingCubes<2>::Edge_Has_Vertex(real phi1,real phi2,real iso_value/*=(real)0*/)
{return (phi1<iso_value&&phi2>=iso_value)||(phi2<iso_value&&phi1>=iso_value);}

MarchingCubes<3>::MarchingCubes(LevelSet<3>& _levelset,std::shared_ptr<SurfaceMesh<3> > _mesh/*=nullptr*/,const real _contour_value/*=(real)0*/)
:levelset(_levelset),grid(_levelset.grid),contour_value(_contour_value)
{
	if(_mesh==nullptr)mesh=std::make_shared<SurfaceMesh<3> >();else mesh=_mesh;

	VectorDi edge_grid_counts[3];
	edge_grid_counts[0]=VectorDi(grid.cell_counts.x()-1,grid.cell_counts.y(),grid.cell_counts.z());
	edge_grid_counts[1]=VectorDi(grid.cell_counts.x(),grid.cell_counts.y()-1,grid.cell_counts.z());
	edge_grid_counts[2]=VectorDi(grid.cell_counts.x(),grid.cell_counts.y(),grid.cell_counts.z()-1);
	for(int i=0;i<3;i++){
		edge_grid[i].Initialize(edge_grid_counts[i],grid.dx);
		edge_grid_array[i].Resize(edge_grid[i].cell_counts);edge_grid_array[i].Fill(-1);}
}

void MarchingCubes<3>::Marching()
{
	const VectorDi& cell_counts=grid.cell_counts-VectorDi::Ones();
	Array<VectorD>& vertices=(*mesh->vertices);vertices.clear();
	Array<VectorDi>& triangles=mesh->elements;triangles.clear();
	for(int i=0;i<cell_counts[0];i++){
		for(int j=0;j<cell_counts[1];j++){
			for(int k=0;k<cell_counts[2];k++){
				VectorDi cell_index(i,j,k);
				unsigned int cell_type=Get_Cell_Type(levelset.phi(VectorDi(i,j,k)),levelset.phi(VectorDi(i+1,j,k)),levelset.phi(VectorDi(i+1,j,k+1)),levelset.phi(VectorDi(i,j,k+1)),
														levelset.phi(VectorDi(i,j+1,k)),levelset.phi(VectorDi(i+1,j+1,k)),levelset.phi(VectorDi(i+1,j+1,k+1)),levelset.phi(VectorDi(i,j+1,k+1)),contour_value);
				////lookup edge table
				int edge_type=edge_table[cell_type];
				for(int ei=0;ei<12;ei++){
					if(Has_Edge(edge_type,ei)&&!Has_Particle_On_Edge(cell_index,ei)){
						////add vertex
						VectorDi v1,v2;Get_Edge_Vertex_Index(cell_index,ei,v1,v2);
						real alpha=(contour_value-levelset.phi(VectorDi(v1[0],v1[1],v1[2])))/(levelset.phi(VectorDi(v2[0],v2[1],v2[2]))-levelset.phi(VectorDi(v1[0],v1[1],v1[2])));
						VectorD pos=(1-alpha)*grid.Center(VectorDi(v1[0],v1[1],v1[2]))+alpha*grid.Center(VectorDi(v2[0],v2[1],v2[2]));
						vertices.push_back(pos);Set_Particle_Index_On_Edge(cell_index,ei,(int)vertices.size()-1);}}
				////lookup triangle table
				for(int ti=0;triangle_table[cell_type][ti]!=-1;ti+=3){
					VectorDi t;
					t[0]=Get_Particle_Index_On_Edge(cell_index,triangle_table[cell_type][ti]);
					t[1]=Get_Particle_Index_On_Edge(cell_index,triangle_table[cell_type][ti+1]);
					t[2]=Get_Particle_Index_On_Edge(cell_index,triangle_table[cell_type][ti+2]);
					triangles.push_back(t);}}}}
	////calculate normals
	if(mesh->normals!=nullptr){
		mesh->normals->resize((int)vertices.size());
		for(auto i=0;i<vertices.size();i++) (*mesh->normals)[i]=levelset.Normal(vertices[i]);}
}

bool MarchingCubes<3>::Has_Edge(const int edge_type,const int edge_index)
{return (edge_type&(1<<edge_index))==0?false:true;}

bool MarchingCubes<3>::Has_Particle_On_Edge(const VectorDi& cell_index,const int edge_index)
{return Get_Particle_Index_On_Edge(cell_index,edge_index)!=-1;}

void MarchingCubes<3>::Get_Edge_Vertex_Index(const VectorDi& cell_index,int edge_index,/*rst*/VectorDi& v1,/*rst*/VectorDi& v2)
{
	switch(edge_index){
		case 0:v1=VectorDi(cell_index[0],cell_index[1],cell_index[2]);v2=VectorDi(cell_index[0]+1,cell_index[1],cell_index[2]);break;
		case 1:v1=VectorDi(cell_index[0]+1,cell_index[1],cell_index[2]);v2=VectorDi(cell_index[0]+1,cell_index[1],cell_index[2]+1);break;
		case 2:v1=VectorDi(cell_index[0]+1,cell_index[1],cell_index[2]+1);v2=VectorDi(cell_index[0],cell_index[1],cell_index[2]+1);break;
		case 3:v1=VectorDi(cell_index[0],cell_index[1],cell_index[2]+1);v2=VectorDi(cell_index[0],cell_index[1],cell_index[2]);break;
		case 4:v1=VectorDi(cell_index[0],cell_index[1]+1,cell_index[2]);v2=VectorDi(cell_index[0]+1,cell_index[1]+1,cell_index[2]);break;
		case 5:v1=VectorDi(cell_index[0]+1,cell_index[1]+1,cell_index[2]);v2=VectorDi(cell_index[0]+1,cell_index[1]+1,cell_index[2]+1);break;
		case 6:v1=VectorDi(cell_index[0]+1,cell_index[1]+1,cell_index[2]+1);v2=VectorDi(cell_index[0],cell_index[1]+1,cell_index[2]+1);break;
		case 7:v1=VectorDi(cell_index[0],cell_index[1]+1,cell_index[2]+1);v2=VectorDi(cell_index[0],cell_index[1]+1,cell_index[2]);break;
		case 8:v1=VectorDi(cell_index[0],cell_index[1],cell_index[2]);v2=VectorDi(cell_index[0],cell_index[1]+1,cell_index[2]);break;
		case 9:v1=VectorDi(cell_index[0]+1,cell_index[1],cell_index[2]);v2=VectorDi(cell_index[0]+1,cell_index[1]+1,cell_index[2]);break;
		case 10:v1=VectorDi(cell_index[0]+1,cell_index[1],cell_index[2]+1);v2=VectorDi(cell_index[0]+1,cell_index[1]+1,cell_index[2]+1);break;
		case 11:v1=VectorDi(cell_index[0],cell_index[1],cell_index[2]+1);v2=VectorDi(cell_index[0],cell_index[1]+1,cell_index[2]+1);break;
		default:std::cerr<<"invalid edge index."<<std::endl;break;}
}

int MarchingCubes<3>::Get_Particle_Index_On_Edge(const VectorDi& cell_index,const int edge_index)
{
	switch(edge_index){
		////axis-x edges
		case 0:return edge_grid_array[0](VectorDi(cell_index[0],cell_index[1],cell_index[2]));
		case 2:return edge_grid_array[0](VectorDi(cell_index[0],cell_index[1],cell_index[2]+1));
		case 4:return edge_grid_array[0](VectorDi(cell_index[0],cell_index[1]+1,cell_index[2]));
		case 6:return edge_grid_array[0](VectorDi(cell_index[0],cell_index[1]+1,cell_index[2]+1));
		////axis-y edges
		case 8:return edge_grid_array[1](VectorDi(cell_index[0],cell_index[1],cell_index[2]));
		case 9:return edge_grid_array[1](VectorDi(cell_index[0]+1,cell_index[1],cell_index[2]));
		case 11:return edge_grid_array[1](VectorDi(cell_index[0],cell_index[1],cell_index[2]+1));
		case 10:return edge_grid_array[1](VectorDi(cell_index[0]+1,cell_index[1],cell_index[2]+1));
		////axis-z edges
		case 3:return edge_grid_array[2](VectorDi(cell_index[0],cell_index[1],cell_index[2]));
		case 1:return edge_grid_array[2](VectorDi(cell_index[0]+1,cell_index[1],cell_index[2]));
		case 7:return edge_grid_array[2](VectorDi(cell_index[0],cell_index[1]+1,cell_index[2]));
		case 5:return edge_grid_array[2](VectorDi(cell_index[0]+1,cell_index[1]+1,cell_index[2]));
		default:std::cerr<<"invalid edge index."<<std::endl;return -1;}
}

void MarchingCubes<3>::Set_Particle_Index_On_Edge(const VectorDi& cell_index,const int edge_index,const int value_input)
{
	switch(edge_index){
		////axis-x edges
		case 0:edge_grid_array[0](VectorDi(cell_index[0],cell_index[1],cell_index[2]))=value_input;break;
		case 2:edge_grid_array[0](VectorDi(cell_index[0],cell_index[1],cell_index[2]+1))=value_input;break;
		case 4:edge_grid_array[0](VectorDi(cell_index[0],cell_index[1]+1,cell_index[2]))=value_input;break;
		case 6:edge_grid_array[0](VectorDi(cell_index[0],cell_index[1]+1,cell_index[2]+1))=value_input;break;
		////axis-y edges
		case 8:edge_grid_array[1](VectorDi(cell_index[0],cell_index[1],cell_index[2]))=value_input;break;
		case 9:edge_grid_array[1](VectorDi(cell_index[0]+1,cell_index[1],cell_index[2]))=value_input;break;
		case 11:edge_grid_array[1](VectorDi(cell_index[0],cell_index[1],cell_index[2]+1))=value_input;break;
		case 10:edge_grid_array[1](VectorDi(cell_index[0]+1,cell_index[1],cell_index[2]+1))=value_input;break;
		////axis-z edges
		case 3:edge_grid_array[2](VectorDi(cell_index[0],cell_index[1],cell_index[2]))=value_input;break;
		case 1:edge_grid_array[2](VectorDi(cell_index[0]+1,cell_index[1],cell_index[2]))=value_input;break;
		case 7:edge_grid_array[2](VectorDi(cell_index[0],cell_index[1]+1,cell_index[2]))=value_input;break;
		case 5:edge_grid_array[2](VectorDi(cell_index[0]+1,cell_index[1]+1,cell_index[2]))=value_input;break;
		default:std::cerr<<"invalid edge index."<<std::endl;break;}
}

unsigned int MarchingCubes<3>::Get_Cell_Type(real v0,real v1,real v2,real v3,real v4,real v5,real v6,real v7,real iso_value)
{
	unsigned int type=0;
	if(v0<iso_value) type|=1;
	if(v1<iso_value) type|=2;
	if(v2<iso_value) type|=4;
	if(v3<iso_value) type|=8;
	if(v4<iso_value) type|=16;
	if(v5<iso_value) type|=32;
	if(v6<iso_value) type|=64;
	if(v7<iso_value) type|=128;
	return type;
}

const int MarchingCubes<3>::edge_table[256] = {
0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

const int MarchingCubes<3>::triangle_table[256][16] =
{{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
{3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
{3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
{3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
{9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
{2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
{8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
{4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
{3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
{1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
{4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
{4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
{5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
{2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
{9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
{0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
{2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
{10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
{5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
{5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
{9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
{1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
{10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
{8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
{2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
{7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
{2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
{11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
{5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
{11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
{11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
{9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
{2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
{6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
{3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
{6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
{10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
{6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
{8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
{7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
{3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
{5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
{0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
{9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
{8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
{5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
{0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
{6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
{10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
{10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
{8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
{1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
{0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
{10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
{3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
{6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
{9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
{8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
{3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
{6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
{0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
{10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
{10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
{2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
{7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
{7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
{2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
{1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
{11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
{8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
{0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
{7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
{10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
{2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
{6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
{7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
{2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
{1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
{10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
{10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
{0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
{7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
{6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
{8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
{9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
{6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
{4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
{10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
{8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
{0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
{1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
{8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
{10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
{4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
{10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
{5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
{11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
{9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
{6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
{7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
{3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
{7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
{3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
{6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
{9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
{1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
{4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
{7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
{6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
{3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
{0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
{6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
{0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
{11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
{6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
{5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
{9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
{1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
{1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
{10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
{0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
{5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
{10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
{11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
{9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
{7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
{2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
{8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
{9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
{9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
{1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
{9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
{9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
{5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
{0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
{10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
{2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
{0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
{0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
{9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
{5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
{3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
{5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
{8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
{0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
{9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
{0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
{1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
{3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
{4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
{9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
{11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
{11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
{2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
{9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
{3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
{1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
{4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
{4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
{0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
{3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
{3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
{0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
{9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
{1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};