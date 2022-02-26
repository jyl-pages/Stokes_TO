#include "EdgeField.h"
#include "File.h"

template<typename T, int d>
inline void EdgeField<T, d>::Resize(const VectorDi & cell_counts)
{
	if (grid.cell_counts == cell_counts) return;
	grid.Initialize(cell_counts);
	for (int i = 0; i < d; i++)
	{
		VectorDi edge_counts = cell_counts + VectorDi::Ones() - VectorDi::Unit(d);
		edge_fields[i].Resize(edge_counts);
	}
}

template<typename T, int d>
inline void EdgeField<T, d>::Fill(const T & value)
{
	for (int i = 0; i < d; i++)
		edge_fields[i].Fill(value);
}

template<typename T, int d>
inline T & EdgeField<T, d>::operator()(const int axis, const VectorDi & coord)
{
	return edge_fields[axis](coord);
}

template<typename T, int d>
inline const T & EdgeField<T, d>::operator()(const int axis, const VectorDi & coord) const
{
	return edge_fields[axis](coord);
}

template<typename T, int d>
inline bool EdgeField<T, d>::Valid_Edge(const int axis, const VectorDi & coord) const
{
	for (int dd = 0; dd < d; dd++)
		if (coord[dd] >= (grid.cell_counts[dd] + (dd != axis)) || coord[dd] < 0) return false;
	return true;
}

template<typename T, int d>
inline void EdgeField<T, d>::Write_Binary(std::ostream & output) const
{
	File::Write_Binary(output, grid.cell_counts);
	for(int i=0;i<d;i++) File::Write_Binary_Array(output, &edge_fields[i].array[0], (int)edge_fields[i].array.size());
}

template<typename T, int d>
inline void EdgeField<T, d>::Write_Binary(const std::string & file_name) const
{
	std::ofstream output(file_name, std::ios::binary);
	if (!output) { std::cerr << "FaceField<T, d>::Write_Binary error: cannot open file " << file_name << "\n "; exit(0); }
	Write_Binary(output);
	output.close();
}

template<typename T, int d>
inline void EdgeField<T, d>::Read_Binary(std::istream & input)
{
	VectorDi cell_counts; File::Read_Binary(input, cell_counts); Resize(cell_counts);
	for (int i = 0; i < d; i++)File::Read_Binary_Array(input, &edge_fields[i].array[0], (int)edge_fields[i].array.size());
}

template<typename T, int d>
inline void EdgeField<T, d>::Read_Binary(const std::string & file_name)
{
	std::ifstream input(file_name, std::ios::binary);
	if (!input) { std::cerr << "FaceField<T, d>::Read_Binary error: cannot open file " << file_name << "\n "; exit(0); }
	Read_Binary(input);
	input.close();

}

template<typename T, int d>
inline void EdgeField<T, d>::Write_To_File_3d(const std::string & file_name, bool extrapolation) const
{
	if constexpr (d == 3) {
		Field<Vector<T, d>, d> v; 
		//Edge_To_Cell_Conversion(*this, v);
		Edge_To_Node_Conversion(*this, v, extrapolation);
		File::Write_Binary_To_File(file_name, v);
	}
	else {
		Field<Vector<T, d>, d> v; 
		//Edge_To_Cell_Conversion(*this, v);
		Edge_To_Node_Conversion(*this, v, extrapolation);
		Field<Vector<T, 3>, 3> v3; 
		VF_Dim_Conversion<T, d, 3>(v, v3);
		File::Write_Binary_To_File(file_name, v3);
	}
}

template<typename T, int d>
void EdgeField<T, d>::Edge_To_Cell_Conversion(const EdgeField<T, d>& edge_field, Field<Vector<T, d>, d>& cell_field)
{
	const Grid<d>& grid = edge_field.grid;
	cell_field.Resize(grid.cell_counts);
	if constexpr (d == 2)
	{
		for (int j = 0; j < grid.cell_counts[1]; j++) for (int i = 0; i < grid.cell_counts[0]; i++)
		{
			T vx(0), vy(0);
			vx = (T)0.5*(edge_field.edge_fields[0](Vector2i(i, j)) + edge_field.edge_fields[0](Vector2i(i, j + 1)));
			vy = (T)0.5*(edge_field.edge_fields[1](Vector2i(i, j)) + edge_field.edge_fields[1](Vector2i(i + 1, j)));
			cell_field(i, j) = Vector<T,d>(vx, vy);
		}
	}
	else
	{
		for (int k = 0; k < grid.cell_counts[2]; k++) for (int j = 0; j < grid.cell_counts[1]; j++) for (int i = 0; i < grid.cell_counts[0]; i++)
		{
			T vx(0), vy(0), vz(0);
			vx = (T)0.25*(
				edge_field.edge_fields[0](Vector3i(i, j, k)) + edge_field.edge_fields[0](Vector3i(i, j + 1, k)) + \
				edge_field.edge_fields[0](Vector3i(i, j, k + 1)) + edge_field.edge_fields[0](Vector3i(i, j + 1, k + 1))
				);
			vy = (T)0.25*(
				edge_field.edge_fields[1](Vector3i(i, j, k)) + edge_field.edge_fields[1](Vector3i(i + 1, j, k)) + \
				edge_field.edge_fields[1](Vector3i(i, j, k + 1)) + edge_field.edge_fields[1](Vector3i(i + 1, j, k + 1))
				);
			vz = (T)0.25*(
				edge_field.edge_fields[2](Vector3i(i, j, k)) + edge_field.edge_fields[2](Vector3i(i + 1, j, k)) + \
				edge_field.edge_fields[2](Vector3i(i, j + 1, k)) + edge_field.edge_fields[2](Vector3i(i + 1, j + 1, k))
				);
			cell_field(i, j, k) = Vector<T, d>(vx, vy, vz);
		}
	}
}

template<typename T, int d>
void EdgeField<T, d>::Edge_To_Node_Conversion(const EdgeField<T, d>& edge_field, Field<Vector<T, d>, d>& node_field, bool extrapolate)
{
	const Grid<d>& grid = edge_field.grid;
	node_field.Resize(grid.node_counts, Vector<T, d>::Zero());
	if constexpr (d == 2)
	{
		int l_x = extrapolate ? 0 : 1, u_x = extrapolate ? grid.node_counts[0] : grid.node_counts[0] - 1;
		int l_y = extrapolate ? 0 : 1, u_y = extrapolate ? grid.node_counts[1] : grid.node_counts[1] - 1;
		for (int j = l_y; j < u_y; j++) for (int i = l_x; i < u_x; i++)
		{
			T v[d];
			for (int dd = 0; dd < d; dd++)
			{
				v[dd] = (T)0;
				int c = 0;

				if (edge_field.Valid_Edge(dd, Vector2i(i, j)))
				{
					v[dd] += edge_field(dd, Vector2i(i, j));
					c++;
				}

				if (edge_field.Valid_Edge(dd, Vector2i(i, j) - Vector2i::Unit(dd)))
				{
					v[dd] += edge_field(dd, Vector2i(i, j) - Vector2i::Unit(dd));
					c++;
				}

				if (c == 0) printf("wtf");
				v[dd] /= c;
			}
			node_field(i, j) = Vector<T, 2>(v[0], v[1]);
		}
	}
	else
	{
		int l_x = extrapolate ? 0 : 1, u_x = extrapolate ? grid.node_counts[0] : grid.node_counts[0] - 1;
		int l_y = extrapolate ? 0 : 1, u_y = extrapolate ? grid.node_counts[1] : grid.node_counts[1] - 1;
		int l_z = extrapolate ? 0 : 1, u_z = extrapolate ? grid.node_counts[2] : grid.node_counts[2] - 1;
		for (int k = l_z; k < u_z; k++) for (int j = l_y; j < u_y; j++) for (int i = l_x; i < u_x; i++)
		{
			T v[d];
			for (int dd = 0; dd < d; dd++)
			{
				v[dd] = (T)0;
				int c = 0;

				if (edge_field.Valid_Edge(dd, Vector3i(i, j, k)))
				{
					v[dd] += edge_field(dd, Vector3i(i, j, k));
					c++;
				}

				if (edge_field.Valid_Edge(dd, Vector3i(i, j, k) - Vector3i::Unit(dd)))
				{
					v[dd] += edge_field(dd, Vector3i(i, j, k) - Vector3i::Unit(dd));
					c++;
				}

				if (c == 0) printf("wtf");
				v[dd] /= c;
			}
			node_field(i, j, k) = Vector<T, 3>(v[0], v[1], v[2]);
		}
	}
}

template class EdgeField<int, 2>;
template class EdgeField<int, 3>;
template class EdgeField<float, 2>;
template class EdgeField<float, 3>;
template class EdgeField<double, 2>;
template class EdgeField<double, 3>;
