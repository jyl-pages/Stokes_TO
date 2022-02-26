#pragma once
#include <fstream>
#include "Field.h"

// minimal functionality
template<typename T, int d>
class EdgeField
{
	Typedef_VectorDi(d);
public:
	Grid<d> grid;
	ArrayF<Field<T, d>, d> edge_fields;

	void Resize(const VectorDi& cell_counts);
	void Fill(const T& value);

	inline T& operator()(const int axis, const VectorDi& coord);
	inline const T& operator()(const int axis, const VectorDi& coord) const;

	inline bool Valid_Edge(const int axis, const VectorDi& coord) const;

	virtual void Write_Binary(std::ostream& output) const;
	virtual void Write_Binary(const std::string& file_name)const;
	virtual void Read_Binary(std::istream& input);
	virtual void Read_Binary(const std::string& file_name);
	virtual void Write_To_File_3d(const std::string& file_name, bool extrapolation=false) const;

	static void Edge_To_Cell_Conversion(const EdgeField<T, d>& edge_field, Field<Vector<T, d>, d>& cell_field);

	static void Edge_To_Node_Conversion(const EdgeField<T, d>& edge_field, Field<Vector<T, d>, d>& node_field, bool extrapolate = false);
};