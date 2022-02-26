//////////////////////////////////////////////////////////////////////////
// Max Grid
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "MacGrid.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////Initializer
template<int d> MacGrid<d>::MacGrid(const VectorDi& _cell_counts,const real _dx,const VectorD& _domain_min)
:grid(_cell_counts,_dx,_domain_min){Initialize_Helper();}

template<int d> MacGrid<d>::MacGrid(const Grid<d>& _grid)
:grid(_grid){Initialize_Helper();}

template<int d> MacGrid<d>& MacGrid<d>::operator=(const MacGrid& copy)
{grid=copy.grid;face_counts=copy.face_counts;face_grids=copy.face_grids;return *this;}

template<int d> void MacGrid<d>::Initialize_Helper()
{Initialize_Face_Counts();Initialize_Face_Grids();}

template<int d> void MacGrid<d>::Initialize_Face_Counts()
{face_counts.resize(d);for(int i=0;i<d;i++){face_counts[i]=grid.cell_counts;face_counts[i][i]+=1;}}

template<int d> void MacGrid<d>::Initialize_Face_Grids()
{face_grids.clear();face_grids.reserve(d);
for(int i=0;i<d;i++){VectorD domain_min=grid.domain_min+(VectorD::Ones()-VectorD::Unit(i))* grid.dx* (real).5;face_grids.emplace_back(face_counts[i]-VectorDi::Ones(),grid.dx,domain_min);}
}

template<int d> void MacGrid<d>::Initialize(const VectorDi& _cell_counts,const real _dx,const VectorD& _domain_min)
{grid.Initialize(_cell_counts,_dx,_domain_min);Initialize_Helper();}

template<int d> void MacGrid<d>::Initialize(const Grid<d>& _grid)
{grid.Initialize(_grid.cell_counts,_grid.dx,_grid.domain_min);Initialize_Helper();}

////////////////////////////////////////////////////////////////////////////////////////////////////
////Face operations
template<int d> bool MacGrid<d>::Valid_Face(const int axis,const VectorDi& face) const
{return face_grids[axis].Valid_Node(face);}

template<int d> bool MacGrid<d>::Is_Boundary_Face(const int axis,const VectorDi& face) const
{return face_grids[axis].Is_Boundary_Node(face);}

template<int d> bool MacGrid<d>::Is_Boundary_Face(const int axis,const VectorDi& face,const int ring_num) const
{
	VectorDi c0=Face_Left_Cell(axis,face);VectorDi c1=Face_Right_Cell(axis,face);
	return ((!grid.Valid_Cell(c0)||grid.Is_Boundary_Cell(c0,ring_num))||(!grid.Valid_Cell(c1)||grid.Is_Boundary_Cell(c1,ring_num)));
}

template<int d> bool MacGrid<d>::Is_Axial_Boundary_Face(const int axis,const VectorDi& face) const
{return face[axis]==0||face[axis]==face_grids[axis].node_counts[axis]-1;}

template<int d> Vector<real,d> MacGrid<d>::Face_Center(const int axis,const VectorDi& face) const
{return face_grids[axis].Node(face);}

template<int d> int MacGrid<d>::Face_Index(const int axis,const VectorDi& face) const
{return face_grids[axis].Node_Index(face);}

template<int d> Vector<int,d> MacGrid<d>::Face_Coord(const int axis,const int face_index) const
{return face_grids[axis].Node_Coord(face_index);}

template<int d> int MacGrid<d>::Number_Of_Faces(const int axis) const
{return face_grids[axis].node_counts.prod();}

template<int d> int MacGrid<d>::Number_Of_Faces() const
{int n=0;for(int axis=0;axis<d;axis++)n+=face_grids[axis].node_counts.prod();return n;}

template<int d> bool MacGrid<d>::Is_Face_Incident_To_Cell(const int axis,const VectorDi& face,const VectorDi& cell)
{return Face_Left_Cell(axis,face)==cell||Face_Right_Cell(axis,face)==cell;}

template<int d> Face<d> MacGrid<d>::Node_Incident_Face(const VectorDi& node,const int i){/*not impl*/return Face<d>();}
template<> Face<2> MacGrid<2>::Node_Incident_Face(const Vector2i& node,const int i)
{static const int axies[4]={0,0,1,1};static const Vector2i face_offset[4]={Vector2i(0,-1),Vector2i(0,0),Vector2i(-1,0),Vector2i(0,0)};
assert(i>=0&&i<4);return Face<2>(axies[i],node+face_offset[i]);}
template<> Face<3> MacGrid<3>::Node_Incident_Face(const Vector3i& node,const int i)
{static const int axies[12]={0,0,0,0,1,1,1,1,2,2,2,2};
static const Vector3i face_offset[12]={Vector3i(0,-1,-1),Vector3i(0,-1,0),Vector3i(0,0,-1),Vector3i(0,0,0),
Vector3i(-1,0,-1),Vector3i(-1,0,0),Vector3i(0,0,-1),Vector3i(0,0,0),Vector3i(-1,-1,0),Vector3i(-1,0,0),Vector3i(0,-1,0),Vector3i(0,0,0)};
assert(i>=0&&i<12);return Face<3>(axies[i],node+face_offset[i]);}

template<int d> Vector<int,d> MacGrid<d>::Node_Incident_Face_Per_Axis(const VectorDi& node,const int i,const int axis)
{int ii=axis*Number_Of_Node_Incident_Faces_Per_Axis()+i;return Node_Incident_Face(node,ii).face;}

template<int d> Vector<int,d> MacGrid<d>::Face_Incident_Node(const int axis,const VectorDi& face,const int i){/*not impl*/return VectorDi::Zero();}
template<> Vector<int,2> MacGrid<2>::Face_Incident_Node(const int axis,const Vector2i& face,const int i)
{static const Vector2i node_offset[4]={Vector2i(0,0),Vector2i(0,1),Vector2i(0,0),Vector2i(1,0)};
assert(i>=0&&i<2);return face+node_offset[axis*2+i];}
template<> Vector<int,3> MacGrid<3>::Face_Incident_Node(const int axis,const Vector3i& face,const int i)
{static const Vector3i node_offset[12]={Vector3i(0,0,0),Vector3i(0,0,1),Vector3i(0,1,0),Vector3i(0,1,1),
Vector3i(0,0,0),Vector3i(0,0,1),Vector3i(1,0,0),Vector3i(1,0,1),Vector3i(0,0,0),Vector3i(0,1,0),Vector3i(1,0,0),Vector3i(1,1,0)};
assert(i>=0&&i<4);return face+node_offset[axis*4+i];}

template<int d> void MacGrid<d>::Face_Incident_Nodes(const int axis,const VectorDi& face,ArrayF2P<VectorDi,d-1>& face_nodes)
{for(int i=0;i<Number_Of_Face_Incident_Nodes();i++)face_nodes[i]=Face_Incident_Node(axis,face,i);}

template<int d> Vector<int,d> MacGrid<d>::Face_Incident_Cell(const int axis,const VectorDi& face,const int i)
{assert(i>=0&&i<2);return face+VectorDi::Unit(axis)*(-1+i);}

template<int d> Vector<int,d> MacGrid<d>::Cell_Incident_Face(const int axis,const VectorDi& cell,const int i)
{assert(i>=0&&i<2);return cell+VectorDi::Unit(axis)*i;}

template<int d> void MacGrid<d>::Cell_Incident_Face(const VectorDi& cell,const int i,int& axis,VectorDi& face)
{
	int side;Grid<d>::Nb_C_Axis_And_Side(i,axis,side);
	face=cell+VectorDi::Unit(axis)*side;	
}

////No instantiation for MacGrid<1>
template class MacGrid<2>;
template class MacGrid<3>;