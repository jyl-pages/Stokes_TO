//////////////////////////////////////////////////////////////////////////
// Grid
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Grid_h__
#define __Grid_h__
#include <fstream>
#include "Common.h"

////Data structure for a uniform Cartersian grid
template<int d> class Grid
{	Typedef_VectorDii(d);
public:
	VectorDi cell_counts;
	VectorDi node_counts;
	real dx;
	VectorD domain_min;
	VectorD domain_max;

	////Initializer
	Grid(const VectorDi& _cell_counts=VectorDi::Zero(),const real _dx=(real)0,const VectorD& _domain_min=VectorD::Zero());
	Grid& operator=(const Grid& copy);
	Grid(const Grid& copy){*this=copy;}
    void Initialize(const VectorDi& _cell_counts=VectorDi::Zero(),const real _dx=(real)0,const VectorD& _domain_min=VectorD::Zero());

	////Coordinate operations
    static inline VectorDi Coord(const int index,const VectorDi& counts);               ////general index->coord
    inline VectorDi Node_Coord(const int index) const;			                        ////node index->node coord
	inline VectorDi Cell_Coord(const int index) const;                                  ////cell index->cell coord
	VectorDi Cell_Coord(const VectorD& pos) const;                                      ////pos->cell coord
	VectorDi Cell_Coord_With_Fraction(const VectorD& pos,VectorD& frac) const;          ////pos->(cell coord,frac), this function has clamp
	VectorDi Cell_Coord_With_Fraction_No_Clamp(const VectorD& pos,VectorD& frac) const;
    
	////Position operations
    VectorD Node(const VectorDi& node) const;                                           ////node coord->node pos
    VectorD Center(const VectorDi& cell) const;                                         ////cell coord->cell center pos
	VectorD Pos(const VectorDi& cell,const VectorD& frac) const;						////cell coord + frac -> interpolated pos                                

    ////Index operations
    static int Index(const Vector<int,d>& coord,const Vector<int,d>& counts);			////general coord->index
    int Cell_Index(const VectorDi& cell) const;                                         ////cell coord->cell index
    int Node_Index(const VectorDi& node) const;                                         ////node coord->node index

    ////Clamp
    VectorD Clamped_Position(const VectorD& pos) const;                                 ////pos->clamped pos
	VectorDi Clamped_Cell_Coord(const VectorD& pos) const;                              ////pos->clamped pos->clamped cell coord
	VectorDi Clamped_Cell_Coord(const VectorDi& cell) const;                            ////cell coord->clamped cell coord
	VectorDi Clamped_Cell_Coord_With_Fraction(const VectorD& pos,VectorD& frac) const;  ////pos->clamped pos->(clamped cell coord,frac)
	VectorDi Clamped_Node_Coord(const VectorDi& node) const;                            ////node coord->clamped node coord

    ////Validate index and coord
	static bool Valid(const VectorDi& coord,const VectorDi& counts);                    ////check if coord is within [zero,counts)
	inline bool Valid_Node(const VectorDi& node) const;                                 ////check if node is within [zero,node_counts)
    inline bool Valid_Cell(const VectorDi& cell) const;									////check if node is within [zero,cell_counts)

	////Counts and sizes
	inline int Number_Of_Cells() const {return cell_counts.prod();}                     ////total number of cells
	inline int Number_Of_Nodes() const {return node_counts.prod();}                     ////total number of nodes
	inline VectorD Length() const {return domain_max-domain_min;}                       ////grid length
	inline VectorD Center() const {return (real).5*(domain_min+domain_max);}            ////grid center
	inline VectorD Position(const VectorD& frac) const {return domain_min+Length().cwiseProduct(frac);}	////relative position inside grid
	VectorD Cell_Frac(const VectorDi& cell) const;                                      ////cell->cell/(cell_counts-1)
	VectorD Node_Frac(const VectorDi& node) const;                                      ////node->node/(node_counts-1)
	bool Inside(const VectorD& pos) const;                                              ////check if pos is inside the grid bounding box with closed boundary

	////Coarsening, refining, and enlarging
	Grid<d> Coarsened(const int factor=2) const;                                        ////return a coarsened grid
    Grid<d> Refined(const int factor=2) const;                                          ////return a refined grid
    Grid<d> Enlarged(const VectorDi& ghost_cell_num_vec) const;                         ////return an enlarged grid
    Grid<d> Enlarged(const int ghost_cell_num) const;                                   ////return an enlarged grid

	////Neighbor operations
    static constexpr int Number_Of_Nb_R(){return Pow(3,d);}								////number of one ring (R) neighbors, 9 in 2D, 27 in 3D
	static int Number_Of_Nb_RV(const VectorDi& coord,const VectorDi& counts);			////number of one ring (R) neighbors valid
	static VectorDi Nb_R(const VectorDi& coord,const int i);							////return coord of R neighbors by index, including itself
	static void Nb_NR(const int n,Array<VectorDi>& cells);								////n ring neighbors, including itself

	static constexpr int Number_Of_Nb_C(){return 2*d;}									////number of one cross (C) neighbors, 4 in 2D, 6 in 3D
	static VectorDi Nb_C(const VectorDi& coord,const int i);							////return coord of C neighbors by index, not including itself
	static int Nb_C_Axis(const int i);													////return C neighbor axis by index
	static void Nb_C_Axis_And_Side(const int i,int& axis,int& side);					////return C neighbor axis and side by index
	
	int Cell_Nb_C_Index(const int cell_idx,const int i) const;							////return a cell's C neighbor index, ONLY WORKS FOR NON-BOUNDARY CELLS!
	int Node_Nb_C_Index(const int node_idx,const int i) const;							////return a node's C neighbor index, ONLY WORKS FOR NON-BOUNDARY NODES!

	static constexpr int Number_Of_Cell_Incident_Nodes(){return Pow(2,d);}				////number of cell incident nodes
	static VectorDi Cell_Incident_Node(const VectorDi& cell_coord,const int i);			////return the coord of a cell incident node by index
	static void Cell_Incident_Nodes(const VectorDi& cell,ArrayF2P<VectorDi,d>& nodes);	////collect the coords of cell incident nodes in a ArrayF2P
	void Cell_Incident_Node_Indices(const VectorDi& cell,ArrayF2P<int,d>& node_indices) const;	////collect indices of cell incident nodes in a ArrayF2P

	static constexpr int Number_Of_Node_Incident_Cells(){return Pow(2,d);}				////number of node incident cells
	static VectorDi Node_Incident_Cell(const VectorDi& node_coord,const int i);			////return the cord of a node incident cell by index

	////Periodic operations
	VectorDi Periodic_Node(const VectorDi& node) const;                                 ////return periodic node coord, e.g., for node (9,8) on a 8x8 grid, the periodic node coord is (0,8)
	VectorDi Periodic_Cell(const VectorDi& cell) const;                                 ////return periodic cell coord
	void Periodic_Cell_Incident_Node_Indices(const VectorDi& cell,ArrayF2P<int,d>& node_indices) const;    ////collect incident periodic node coords of a cell

	////Neighbor index accelerator helper
protected:
	ArrayF<int,d*2> cell_nb_c_idx_offset;
	ArrayF<int,d*2> node_nb_c_idx_offset;
	void Initialize_Nb_C_Index_Offset();

	////Boundary operations
protected:
	static bool Is_Boundary(const VectorDi& coord,const VectorDi& domain_min,const VectorDi& domain_max);
	static bool Is_Boundary(const VectorDi& coord,const VectorDi& counts);
	static bool Is_Boundary(const VectorDi& coord,const VectorDi& counts,const int ring_num);
	static bool Is_Boundary_Side(const VectorDi& coord,const VectorDi& domain_min,const VectorDi& domain_max,int& side);
	static bool Is_Axial_Boundary(const VectorDi& coord,const VectorDi& counts,const int ring_num,const int axis);
	static bool Is_Axial_Boundary_Side(const int axis,const VectorDi& coord,const VectorDi& domain_min,const VectorDi& domain_max,int& side);
	static bool Is_Corner(const VectorDi& coord,const VectorDi& counts);
	static bool Is_Corner(const VectorDi& coord,const VectorDi& domain_min,const VectorDi& domain_max);
public:
	bool Is_Boundary_Node(const VectorDi& node) const {return Is_Boundary(node,node_counts);}                               ////check if a node is on the grid boundary
	bool Is_Boundary_Node(const VectorDi& node,const int ring_num) const {return Is_Boundary(node,node_counts,ring_num);}   ////check if a node is on the n-ring grid boundary
	bool Is_Boundary_Cell(const VectorDi& cell) const {return Is_Boundary(cell,cell_counts);}                               ////check if a cell is on the grid boundary
	bool Is_Boundary_Cell(const VectorDi& cell,const int ring_num) const {return Is_Boundary(cell,cell_counts,ring_num);}   ////check if a cell is on the n-ring grid boundary
	bool Is_Corner_Node(const VectorDi& node) const {return Is_Corner(node,node_counts);}	////check if a node is on the grid corner

	////IO
	void Write_Binary(std::ostream& output) const;
	void Read_Binary(std::istream& input);
	void Write_To_File_3d(const std::string& file_name) const;

	/////Conversion
	Grid<d> Node_Grid_To_Cell_Grid() const {return Grid<d>(cell_counts+VectorDi::Ones(),dx,domain_min-dx*VectorD::Ones()*(real).5);}
	Grid<d> Cell_Grid_To_Node_Grid() const {return Grid<d>(cell_counts-VectorDi::Ones(),dx,domain_min+dx*VectorD::Ones()*(real).5);}
};

//////////////////////////////////////////////////////////////////////////
////inline member functions
////Coord operations
template<> inline Vector<int,1> Grid<1>::Coord(const int index,const Vector<int,1>& counts) 
{Vector1i coord;coord[0]=index;return coord;}

template<> inline Vector<int,2> Grid<2>::Coord(const int index,const Vector<int,2>& counts) 
{int x=index/counts[1];int y=index-x*counts[1];return Vector2i(x,y);}

template<> inline Vector<int,3> Grid<3>::Coord(const int index,const Vector<int,3>& counts) 
{
	//int n_yz=counts[1]*counts[2];int mod_yz=index%n_yz;
	//return Vector3i(index/n_yz,mod_yz/counts[2],mod_yz%counts[2]);
	int n_yz=counts[1]*counts[2];
	int x=index/n_yz;
	int mod_yz=index-x*n_yz;
	int y=mod_yz/counts[2];
	int z=mod_yz-y*counts[2];
	return Vector3i(x,y,z);
}

template<int d> inline Vector<int,d> Grid<d>::Node_Coord(const int index) const {return Coord(index,node_counts);}

template<int d> inline Vector<int,d> Grid<d>::Cell_Coord(const int index) const {return Coord(index,cell_counts);}

////Index operations
template<> inline int Grid<1>::Index(const Vector<int,1>& coord,const Vector<int,1>& counts)
{return coord[0];}

template<> inline int Grid<2>::Index(const Vector<int,2>& coord,const Vector<int,2>& counts)
{return coord[0]*counts[1]+coord[1];}

template<> inline int Grid<3>::Index(const Vector<int,3>& coord,const Vector<int,3>& counts)
{return coord[0]*counts[1]*counts[2]+coord[1]*counts[2]+coord[2];}

template<int d> inline int Grid<d>::Cell_Index(const VectorDi& cell) const {return Index(cell,cell_counts);}

template<int d> inline int Grid<d>::Node_Index(const VectorDi& node) const {return Index(node,node_counts);}

////Validation operations
template<> inline bool Grid<1>::Valid_Node(const Vector1i& node) const 
{return node[0]>=0&&node[0]<node_counts[0];}

template<> inline bool Grid<2>::Valid_Node(const Vector2i& node) const 
{return node[0]>=0&&node[0]<node_counts[0]&&node[1]>=0&&node[1]<node_counts[1];}

template<> inline bool Grid<3>::Valid_Node(const Vector3i& node) const 
{return node[0]>=0&&node[0]<node_counts[0]&&node[1]>=0&&node[1]<node_counts[1]&&node[2]>=0&&node[2]<node_counts[2];}

template<> inline bool Grid<1>::Valid_Cell(const Vector1i& cell) const 
{return cell[0]>=0&&cell[0]<cell_counts[0];}

template<> inline bool Grid<2>::Valid_Cell(const Vector2i& cell) const 
{return cell[0]>=0&&cell[0]<cell_counts[0]&&cell[1]>=0&&cell[1]<cell_counts[1];}

template<> inline bool Grid<3>::Valid_Cell(const Vector3i& cell) const 
{return cell[0]>=0&&cell[0]<cell_counts[0]&&cell[1]>=0&&cell[1]<cell_counts[1]&&cell[2]>=0&&cell[2]<cell_counts[2];}

////Dimension conversion
template<int d1,int d2> void Dim_Conversion(const Grid<d1>& input,Grid<d2>& output);

////Lattice functions
template<int d> bool Valid_Node(const Grid<d>& lattice,const Vector<int,d>& node,const Array<int>& mat_array);

////Iterator for traversing the nodes or cells
template<int d> class GridIterator
{Typedef_VectorDi(d);
public:
    enum class IteratorType : int {Node=0,Cell};
    VectorDi start=VectorDi::Zero();
    VectorDi counts=VectorDi::Zero();
    VectorDi coord=VectorDi::Zero();
    int size=0;
    int index=0;

    GridIterator(const VectorDi& _counts,const VectorDi& _start=VectorDi::Zero());
    GridIterator(const Grid<d>& grid,const IteratorType type=IteratorType::Node);

    bool Valid() const {return index<size;}
    void Next();
    int Index() const {return index;}
    VectorDi Coord() const {return start+coord;}
};

////Macros for using GridIterator
#define iterate_cell(iter,grid) \
for(GridIterator<d> iter(grid,GridIterator<d>::IteratorType::Cell);iter.Valid();iter.Next())
#define iterate_node(iter,grid) \
for(GridIterator<d> iter(grid,GridIterator<d>::IteratorType::Node);iter.Valid();iter.Next())
#define iterate_cell_d(iter,grid,dim) \
for(GridIterator<dim> iter(grid,GridIterator<dim>::IteratorType::Cell);iter.Valid();iter.Next())
#define iterate_node_d(iter,grid,dim) \
for(GridIterator<dim> iter(grid,GridIterator<dim>::IteratorType::Node);iter.Valid();iter.Next())
////iterate a subgrid, *including* the end
#define iterate_range(iter,start,end) \
for(GridIterator<d> iter(end-start+Vector<int,d>::Ones(),start);iter.Valid();iter.Next())
#define iterate_range_d(iter,start,end,dim) \
for(GridIterator<dim> iter(end-start+Vector<int,dim>::Ones(),start);iter.Valid();iter.Next())

////CODESAMPLE: openmp version of iterate_cell
//int cell_num=grid.Number_Of_Cells();
//#pragma omp parallel for
//for(int i=0;i<cell_num;i++){const VectorDi& cell=grid.Cell_Coord(i);

////CODESAMPLE: openmp version of iterate_node
//int node_num=grid.Number_Of_Nodes();
//#pragma omp parallel for
//for(int i=0;i<node_num;i++){const VectorDi& node=grid.Node_Coord(i);
#endif