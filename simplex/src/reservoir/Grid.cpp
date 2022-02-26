//////////////////////////////////////////////////////////////////////////
// Grid
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Grid.h"
#include "AuxFunc.h"
#include "File.h"

using namespace AuxFunc;
////////////////////////////////////////////////////////////////////////////////////////////////////
////Initializer
template<int d> Grid<d>::Grid(const VectorDi& _cell_counts,const real _dx,const VectorD& _domain_min)
:cell_counts(_cell_counts),node_counts(_cell_counts+VectorDi::Ones()),dx(_dx),domain_min(_domain_min)
{domain_max=domain_min+cell_counts.template cast<real>()*dx;Initialize_Nb_C_Index_Offset();}

template<int d> Grid<d>& Grid<d>::operator=(const Grid& copy)
{cell_counts=copy.cell_counts;node_counts=copy.node_counts;dx=copy.dx;domain_min=copy.domain_min;domain_max=copy.domain_max;
cell_nb_c_idx_offset=copy.cell_nb_c_idx_offset;node_nb_c_idx_offset=copy.node_nb_c_idx_offset;
return *this;}

template<int d> void Grid<d>::Initialize(const VectorDi& _cell_counts,const real _dx,const VectorD& _domain_min)
{cell_counts=_cell_counts;node_counts=cell_counts+VectorDi::Ones();dx=_dx;domain_min=_domain_min;domain_max=domain_min+cell_counts.template cast<real>()*dx;
Initialize_Nb_C_Index_Offset();}

template<int d> void Grid<d>::Initialize_Nb_C_Index_Offset()
{
	if constexpr (d==1){
		cell_nb_c_idx_offset={-1,1};
		node_nb_c_idx_offset={-1,1};}
	else if constexpr (d==2){
		int y=cell_counts[1];
		cell_nb_c_idx_offset={-y,-1,1,y};
		int ny=node_counts[1];
		node_nb_c_idx_offset={-ny,-1,1,ny};}
	else if constexpr (d==3){
		int y=cell_counts[1];int z=cell_counts[2];
		cell_nb_c_idx_offset={-y*z,-z,-1,1,z,y*z};
		int ny=node_counts[1];int nz=node_counts[2];
		node_nb_c_idx_offset={-ny*nz,-nz,-1,1,nz,ny*nz};}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////Coordinate operations
template<int d> Vector<int,d> Grid<d>::Cell_Coord(const VectorD& pos) const {VectorD coord_with_frac=(pos-domain_min)/dx;return coord_with_frac.template cast<int>();}

template<int d> Vector<int,d> Grid<d>::Cell_Coord_With_Fraction(const VectorD& pos,VectorD& frac) const
{VectorD coord_with_frac=(pos-domain_min)/dx;VectorDi coord=coord_with_frac.template cast<int>();coord=Clamped_Cell_Coord(coord);frac=coord_with_frac-coord.template cast<real>();return coord;}

template<int d> Vector<int,d> Grid<d>::Cell_Coord_With_Fraction_No_Clamp(const VectorD& pos,VectorD& frac) const
{
	VectorD coord_with_frac=(pos-domain_min)/dx;
	for(int i=0;i<d;i++)frac[i]=std::floor(coord_with_frac[i]);
	VectorDi coord=frac.template cast<int>();
	frac=coord_with_frac-frac;return coord;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////Position operations
template<int d> Vector<real,d> Grid<d>::Node(const VectorDi& node) const {return domain_min+node.template cast<real>()*dx;}

template<int d> Vector<real,d> Grid<d>::Center(const VectorDi& cell) const {return domain_min+(cell.template cast<real>()+(real).5*VectorD::Ones())*dx;}

template<int d> Vector<real,d> Grid<d>::Pos(const VectorDi& cell,const VectorD& frac) const {return Node(cell)+frac*dx;}

////Clamp
template<int d> Vector<real,d> Grid<d>::Clamped_Position(const VectorD& pos) const     ////Attention: clamp with epsilon on the max side
{VectorD domain_max_minus_epsilon=domain_max-VectorD::Ones()*dx*(real)1e-6;return Clamp_Vector(pos,domain_min,domain_max_minus_epsilon);}

template<int d> Vector<int,d> Grid<d>::Clamped_Cell_Coord(const VectorDi& cell) const
{const VectorDi lower=VectorDi::Zero();const VectorDi upper=cell_counts-VectorDi::Ones();return Clamp_Vector(cell,lower,upper);}

template<int d> Vector<int,d> Grid<d>::Clamped_Cell_Coord(const VectorD& pos) const
{return Clamped_Cell_Coord(Cell_Coord(pos));}

template<int d> Vector<int,d> Grid<d>::Clamped_Cell_Coord_With_Fraction(const VectorD& pos,VectorD& frac) const
{VectorD clamped_pos=Clamped_Position(pos);return Cell_Coord_With_Fraction(clamped_pos,frac);}

template<int d> Vector<int,d> Grid<d>::Clamped_Node_Coord(const VectorDi& node) const
{const VectorDi lower=VectorDi::Zero();const VectorDi upper=node_counts-VectorDi::Ones();return Clamp_Vector(node,lower,upper);}

////////////////////////////////////////////////////////////////////////////////////////////////////
////Index operations

////////////////////////////////////////////////////////////////////////////////////////////////////
////Validate index and coord
template<int d> bool Grid<d>::Valid(const VectorDi& coord,const VectorDi& counts)
{static const VectorDi zero=VectorDi::Zero();return All_Greater_Equal(coord,zero)&&All_Less(coord,counts);}

////Definitions of Valid_Node and Valid_Cell are in the header file

////////////////////////////////////////////////////////////////////////////////////////////////////
////Counts and sizes
template<int d> Vector<real,d> Grid<d>::Cell_Frac(const VectorDi& cell) const
{VectorD f;for(int i=0;i<d;i++)f[i]=(real)cell[i]/(real)(cell_counts[i]-1);return f;}

template<int d> Vector<real,d> Grid<d>::Node_Frac(const VectorDi& node) const
{VectorD f;for(int i=0;i<d;i++)f[i]=(real)node[i]/(real)(node_counts[i]-1);return f;}

template<int d> bool Grid<d>::Inside(const VectorD& pos) const
{for(int i=0;i<d;i++)if(pos[i]<domain_min[i]||pos[i]>domain_max[i])return false;return true;}

////////////////////////////////////////////////////////////////////////////////////////////////////
//Coarsening, refining, enlarging
template<int d> Grid<d> Grid<d>::Coarsened(const int factor) const
{
    Grid<d> coarsened_grid=*this;coarsened_grid.cell_counts/=factor;
    bool has_zero_dim=false;for(int i=0;i<d;i++)if(coarsened_grid.cell_counts[i]==0){coarsened_grid.cell_counts[i]=1;has_zero_dim=true;}
    coarsened_grid.node_counts=coarsened_grid.cell_counts+VectorDi::Ones();coarsened_grid.dx*=(real)factor;
    if(has_zero_dim){
        VectorD center=coarsened_grid.Center();
        VectorD coarse_domain_size=coarsened_grid.cell_counts.template cast<real>()*coarsened_grid.dx;
        coarsened_grid.domain_min=center-(real).5*coarse_domain_size;
        coarsened_grid.domain_max=center+(real).5*coarse_domain_size;}
    return coarsened_grid;
}

template<int d> Grid<d> Grid<d>::Refined(const int factor/*=2*/) const
{const VectorDi refined_cell_counts=cell_counts*factor;real refined_dx=dx/(real)factor;
Grid<d> refined_grid(refined_cell_counts,refined_dx,domain_min);return refined_grid;}

template<int d> Grid<d> Grid<d>::Enlarged(const VectorDi& ghost_cell_num_vec) const
{VectorDi ghost_cell_counts=cell_counts+ghost_cell_num_vec*2;
VectorD ghost_domain_min=domain_min-ghost_cell_num_vec.template cast<real>()*dx;
Grid<d> ghost_Grid;ghost_Grid.Initialize(ghost_cell_counts,dx,ghost_domain_min);return ghost_Grid;}

template<int d> Grid<d> Grid<d>::Enlarged(const int ghost_cell_num) const
{VectorDi ghost_cell_num_vec=VectorDi::Ones()*ghost_cell_num;return Enlarged(ghost_cell_num_vec);}

////////////////////////////////////////////////////////////////////////////////////////////////////
////Neighbor operations
int Number_Of_Nb_RV_1D(const int coord,const int count){int n=3;if(coord==0)n--;if(coord==count-1)n--;return n;}

template<> int Grid<1>::Number_Of_Nb_RV(const Vector1i& coord,const Vector1i& counts)
{return Number_Of_Nb_RV_1D(coord[0],counts[0]);}

template<> int Grid<2>::Number_Of_Nb_RV(const Vector2i& coord,const Vector2i& counts)
{return Number_Of_Nb_RV_1D(coord[0],counts[0])*Number_Of_Nb_RV_1D(coord[1],counts[1]);}

template<> int Grid<3>::Number_Of_Nb_RV(const Vector3i& coord,const Vector3i& counts)
{return Number_Of_Nb_RV_1D(coord[0],counts[0])*Number_Of_Nb_RV_1D(coord[1],counts[1])*Number_Of_Nb_RV_1D(coord[2],counts[2]);}

template<int d> Vector<int,d> Grid<d>::Nb_R(const VectorDi& coord,const int i){/*not impl*/return VectorDi::Zero();}
template<> Vector1i Grid<1>::Nb_R(const Vector1i& coord,const int index)
{assert(index>=0&&index<3);return coord+Vector1i(-1+index);}
template<> Vector2i Grid<2>::Nb_R(const Vector2i& coord,const int index)
{assert(index>=0&&index<9);int i=index/3;int j=index%3;return coord+Vector2i(-1+i,-1+j);}
template<> Vector3i Grid<3>::Nb_R(const Vector3i& coord,const int index)
{assert(index>=0&&index<27);int i=index/9;int m=index%9;int j=m/3;int k=m%3;return coord+Vector3i(-1+i,-1+j,-1+k);}

template<int d> void Grid<d>::Nb_NR(const int n,Array<VectorDi>& cells){/*not impl*/}
template<> void Grid<1>::Nb_NR(const int n,Array<Vector1i>& cells)
{cells.clear();for(int i=-n;i<=n;i++)cells.push_back(Vec1i(i));}
template<> void Grid<2>::Nb_NR(const int n,Array<Vector2i>& cells)
{cells.clear();for(int i=-n;i<=n;i++)for(int j=-n;j<=n;j++)cells.push_back(Vector2i(i,j));}
template<> void Grid<3>::Nb_NR(const int n,Array<Vector3i>& cells)
{cells.clear();for(int i=-n;i<=n;i++)for(int j=-n;j<=n;j++)for(int k=-n;k<=n;k++)cells.push_back(Vector3i(i,j,k));}

template<int d> Vector<int,d> Grid<d>::Nb_C(const VectorDi& coord,const int i){/*not impl*/return VectorDi::Zero();}
template<> Vector<int,1> Grid<1>::Nb_C(const Vector1i& coord,const int i)
{static int neighbor_offset[2]={-1,1};return coord+Vec1i(neighbor_offset[i]);}
template<> Vector<int,2> Grid<2>::Nb_C(const Vector2i& coord,const int i)
{static Vector2i neighbor_offset[4]={Vector2i(-1,0),Vector2i(0,-1),Vector2i(0,1),Vector2i(1,0)};assert(i>=0&&i<4);return coord+neighbor_offset[i];}
template<> Vector<int,3> Grid<3>::Nb_C(const Vector3i& coord,const int i)
{static const Vector3i neighbor_offset[6]={Vector3i(-1,0,0),Vector3i(0,-1,0),Vector3i(0,0,-1),Vector3i(0,0,1),Vector3i(0,1,0),Vector3i(1,0,0)};
assert(i>=0&&i<6);return coord+neighbor_offset[i];}

template<int d> int Grid<d>::Nb_C_Axis(const int i){return 0;}
template<> int Grid<2>::Nb_C_Axis(const int i){static int axies[4]={0,1,1,0};return axies[i];}
template<> int Grid<3>::Nb_C_Axis(const int i){static int axies[6]={0,1,2,2,1,0};return axies[i];}

template<int d> void Grid<d>::Nb_C_Axis_And_Side(const int i,int& axis,int& side){/*not impl*/}
template<> void Grid<1>::Nb_C_Axis_And_Side(const int i,int& axis,int& side){static int sides[2]={0,1};axis=0;side=sides[i];}
template<> void Grid<2>::Nb_C_Axis_And_Side(const int i,int& axis,int& side){static int axies[4]={0,1,1,0};static int sides[4]={0,0,1,1};axis=axies[i];side=sides[i];}
template<> void Grid<3>::Nb_C_Axis_And_Side(const int i,int& axis,int& side){static int axies[6]={0,1,2,2,1,0};static int sides[6]={0,0,0,1,1,1};axis=axies[i];side=sides[i];}


template<int d> int Grid<d>::Cell_Nb_C_Index(const int cell_idx,const int i) const {return cell_idx+cell_nb_c_idx_offset[i];}

template<int d> int Grid<d>::Node_Nb_C_Index(const int node_idx, const int i) const {return node_idx+node_nb_c_idx_offset[i];}

template<int d> Vector<int,d> Grid<d>::Cell_Incident_Node(const VectorDi& cell_coord,const int i){return VectorDi::Zero();}
template<> Vector<int,1> Grid<1>::Cell_Incident_Node(const Vector1i& cell_coord,const int index)
{int i=index/2;return cell_coord+Vec1i(i);}

template<> Vector<int,2> Grid<2>::Cell_Incident_Node(const Vector2i& cell_coord,const int index)
{int i=index/2;int j=index%2;return cell_coord+Vector2i(i,j);}

template<> Vector<int,3> Grid<3>::Cell_Incident_Node(const Vector3i& cell_coord,const int index)
{int i=index/4;int j=(index%4)/2;int k=index%2;return cell_coord+Vector3i(i,j,k);}

template<int d> void Grid<d>::Cell_Incident_Nodes(const VectorDi& cell,ArrayF2P<VectorDi,d>& nodes)
{for(int i=0;i<Grid<d>::Number_Of_Cell_Incident_Nodes();i++)nodes[i]=Grid<d>::Cell_Incident_Node(cell,i);}

template<int d> void Grid<d>::Cell_Incident_Node_Indices(const VectorDi& cell,ArrayF2P<int,d>& node_indices) const
{for(int i=0;i<Grid<d>::Number_Of_Cell_Incident_Nodes();i++){const VectorDi& node=Grid<d>::Cell_Incident_Node(cell,i);node_indices[i]=Node_Index(node);}}

template<int d> Vector<int,d> Grid<d>::Node_Incident_Cell(const VectorDi& node_coord,const int i)
{return Cell_Incident_Node(node_coord,i)-VectorDi::Ones();}

////////////////////////////////////////////////////////////////////////////////////////////////////
////Periodic operations
template<int d> Vector<int,d> Grid<d>::Periodic_Node(const VectorDi& node) const
{VectorDi p_node=node;for(int i=0;i<d;i++){if(node[i]<0)p_node[i]+=node_counts[i];else if(node[i]>=node_counts[i])p_node[i]-=node_counts[i];}return p_node;}

template<int d> Vector<int,d> Grid<d>::Periodic_Cell(const VectorDi& cell) const
{VectorDi p_cell=cell;for(int i=0;i<d;i++){if(cell[i]<0)p_cell[i]+=cell_counts[i];else if(cell[i]>=cell_counts[i])p_cell[i]-=cell_counts[i];}return p_cell;}

template<int d> void Grid<d>::Periodic_Cell_Incident_Node_Indices(const VectorDi& cell,ArrayF2P<int, d>& node_indices) const
{for(int i=0;i<Number_Of_Cell_Incident_Nodes();i++){const VectorDi& node=Periodic_Node(Grid<d>::Cell_Incident_Node(cell,i));node_indices[i]=Node_Index(node);}}

////////////////////////////////////////////////////////////////////////////////////////////////////
////Boundary operations
template<int d> bool Grid<d>::Is_Boundary(const VectorDi& coord,const VectorDi& domain_min,const VectorDi& domain_max)
{return Has_Equal(coord,domain_min)||Has_Equal(coord,domain_max);}

template<int d> bool Grid<d>::Is_Boundary(const VectorDi& coord,const VectorDi& counts)
{static VectorDi zero=VectorDi::Zero();VectorDi coord_max=counts-VectorDi::Ones();return Has_Equal(coord,zero)||Has_Equal(coord,coord_max);}

template<int d> bool Grid<d>::Is_Boundary(const VectorDi& coord,const VectorDi& counts,const int ring_num)
{VectorDi zero=VectorDi::Zero()+VectorDi::Ones()*(ring_num-1);VectorDi coord_max=counts-VectorDi::Ones()*ring_num;return Has_Less_Equal(coord,zero)||Has_Greater_Equal(coord,coord_max);}

template<int d> bool Grid<d>::Is_Axial_Boundary(const VectorDi& coord,const VectorDi& counts,const int ring_num,const int axis)
{int zero=ring_num-1;int coord_max=counts[axis]-ring_num;return coord[axis]<=zero||coord[axis]>=coord_max;}

template<int d> bool Grid<d>::Is_Axial_Boundary_Side(const int axis,const VectorDi& coord,const VectorDi& domain_min,const VectorDi& domain_max,int& side)
{bool left=coord[axis]==domain_min[axis];bool right=coord[axis]==domain_max[axis];
if(left)side=0;else if(right)side=1;return left||right;}

template<int d> bool Grid<d>::Is_Boundary_Side(const VectorDi& coord,const VectorDi& domain_min,const VectorDi& domain_max,int& side)
{bool left=Has_Equal(coord,domain_min);bool right=Has_Equal(coord,domain_max);if(left)side=0;else if(right) side=1;return left||right;}

template<int d> bool Grid<d>::Is_Corner(const VectorDi& coord,const VectorDi& counts)
{for(int i=0;i<d;i++)if(!Is_Axial_Boundary(coord,counts,1,i))return false;return true;}

template<int d> bool Grid<d>::Is_Corner(const VectorDi& coord,const VectorDi& domain_min,const VectorDi& domain_max)
{int side=0;for(int i=0;i<d;i++)if(!Is_Axial_Boundary_Side(i,coord,domain_min,domain_max,side))return false;return true;}

////////////////////////////////////////////////////////////////////////////////////////////////////
////IO
template<int d> void Grid<d>::Write_Binary(std::ostream& output) const
{
    File::Write_Binary(output,cell_counts);
    File::Write_Binary(output,dx);
    File::Write_Binary(output,domain_min);
}

template<int d> void Grid<d>::Read_Binary(std::istream& input)
{
    File::Read_Binary(input,cell_counts);
    File::Read_Binary(input,dx);
    File::Read_Binary(input,domain_min);
    Initialize(cell_counts,dx,domain_min);
}

template<int d> void Grid<d>::Write_To_File_3d(const std::string& file_name) const
{
	if constexpr (d==3)
		File::Write_Binary_To_File(file_name,*this);
	else{
		Grid<3> g3;Dim_Conversion<d,3>(*this,g3);
		File::Write_Binary_To_File(file_name,g3);}
}

template class Grid<1>;
template class Grid<2>;
template class Grid<3>;

////Dimension conversion
template<int d1,int d2> void Dim_Conversion(const Grid<d1>& input,Grid<d2>& output)
{
    constexpr int n=d1<d2?d1:d2;
    Vector<int,d2> cell_counts_output=Vector<int,d2>::Zero();for(int i=0;i<n;i++)cell_counts_output[i]=input.cell_counts[i];
    Vector<real,d2> domain_min_output=Vector<real,d2>::Zero();for(int i=0;i<n;i++)domain_min_output[i]=input.domain_min[i];
    output=Grid<d2>(cell_counts_output,input.dx,domain_min_output);
}

template void Dim_Conversion<2,2>(const Grid<2>&,Grid<2>&);
template void Dim_Conversion<3,3>(const Grid<3>&,Grid<3>&);
template void Dim_Conversion<2,3>(const Grid<2>&,Grid<3>&);
template void Dim_Conversion<3,2>(const Grid<3>&,Grid<2>&);

////Lattice functions
template<int d> bool Valid_Node(const Grid<d>& lattice,const Vector<int,d>& node,const Array<int>& mat_array)
{for(int ii=0;ii<Grid<d>::Number_Of_Node_Incident_Cells();ii++){
        Vector<int,d> cell=lattice.Node_Incident_Cell(node,ii);
        if(lattice.Valid_Cell(cell)&&mat_array[lattice.Cell_Index(cell)]!=-1)return true;}return false;}

template bool Valid_Node<2>(const Grid<2>&,const Vector2i&,const Array<int>&);
template bool Valid_Node<3>(const Grid<3>&,const Vector3i&,const Array<int>&);

template<int d> GridIterator<d>::GridIterator(const VectorDi& _counts,const VectorDi& _start)
:start(_start),counts(_counts),size(counts.prod()){}

template<int d> GridIterator<d>::GridIterator(const Grid<d>& grid,const IteratorType type)
{
    switch (type){
    case IteratorType::Node:{counts=grid.node_counts;}break;
    case IteratorType::Cell:{counts=grid.cell_counts;}break;}
    size=counts.prod();
}

template<> void GridIterator<1>::Next()
{
    coord[0]++;index++;
}

template<> void GridIterator<2>::Next()
{
    coord[1]++;index++;
    if(coord[1]>=counts[1]){coord[1]=0;coord[0]++;}
}

template<> void GridIterator<3>::Next()
{
    coord[2]++;index++;
    if(coord[2]>=counts[2]){coord[2]=0;coord[1]++;if(coord[1]>=counts[1]){coord[1]=0;coord[0]++;}}
}

template class GridIterator<1>;
template class GridIterator<2>;
template class GridIterator<3>;