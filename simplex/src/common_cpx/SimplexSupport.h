#pragma once

#include "MacGrid.h"
#include "Grid.h"
#include "Field.h"
#include "EdgeField.h"
#include "FaceField.h"

#include "form01mapping.h"
#include "grid3D.h"

namespace simplexSupport
{
	enum FieldType
	{
		Cell,
		Node
	};

	// 2D

	template<typename T> 
	void simplex2solver(const Field<T, 2>& x, const grid2D& grid_ref, T *y);

	template<typename T>
	void simplex2solver(const FaceField<T, 2>& x, const grid2D& grid_ref, T *y);

	template<typename T>
	void solver2simplex(const grid2D& grid_ref, T *x, Field<T, 2>& y);

	template<typename T>
	void solver2simplex(const grid2D& grid_ref, T *x, FaceField<T, 2>& y);

	template<typename A,typename B,typename F>
	void simplex2solver(const Field<A, 2>& x, const grid2D& grid_ref, B *y, F f);

	template<typename A, typename B, typename F>
	void simplex2solver(const FaceField<A, 2>& x, const grid2D& grid_ref, B *y, F f);

	template<typename A, typename B, typename F>
	void solver2simplex(const grid2D& grid_ref, A *x, Field<B, 2>& y, F f);

	template<typename A, typename B, typename F>
	void solver2simplex(const grid2D& grid_ref, A *x, FaceField<B, 2>& y, F f);

	// 3D

	template<typename T>
	void simplex2solver(const Field<T, 3>& x, const grid3D& grid_ref, T *y, FieldType type = Cell);

	template<typename T>
	void simplex2solver(const EdgeField<T, 3>&x, const grid3D& grid_ref, T *y);

	template<typename T>
	void simplex2solver(const FaceField<T, 3>& x, const grid3D& grid_ref, T *y);

	template<typename T>
	void solver2simplex(const grid3D& grid_ref, T *x, Field<T, 3>& y, FieldType type = Cell);

	template<typename T>
	void solver2simplex(const grid3D& grid_ref, T *x, EdgeField<T, 3>& y);

	template<typename T>
	void solver2simplex(const grid3D& grid_ref, T *x, FaceField<T, 3>& y);

	template<typename A, typename B, typename F>
	void simplex2solver(const Field<A, 3>& x, const grid3D& grid_ref, B *y, F f, FieldType type = Cell);

	template<typename A, typename B, typename F>
	void simplex2solver(const EdgeField<A, 3>& x, const grid3D& grid_ref, B *y, F f);

	template<typename A, typename B, typename F>
	void simplex2solver(const FaceField<A, 3>& x, const grid3D& grid_ref, B *y, F f);

	template<typename A, typename B, typename F>
	void solver2simplex(const grid3D& grid_ref, A *x, Field<B, 3>& y, F f, FieldType type = Cell);

	template<typename A, typename B, typename F>
	void solver2simplex(const grid3D& grid_ref, A *x, EdgeField<B, 3>& y, F f);

	template<typename A, typename B, typename F>
	void solver2simplex(const grid3D& grid_ref, A *x, FaceField<B, 3>& y, F f);
};

template<typename T>
void simplexSupport::simplex2solver(const Field<T, 2>& x, const grid2D & grid_ref, T * y)
{
	const int cell_off = 0;

	memset(y+cell_off, 0, sizeof(T)*grid_ref.cell_size());

	for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
	{
		int cell_ind = cell_off + grid_ref.cell_ind(i, j);
		Vector2i cell = Vector2i(i, j);
		y[cell_ind] = x(cell);
	}
}

template<typename T>
void simplexSupport::simplex2solver(const FaceField<T, 2>& x, const grid2D & grid_ref, T *y)
{
	const int face_x_off = 0;
	const int face_y_off = face_x_off + grid_ref.face_size(0);

	memset(y, 0, sizeof(T)*(grid_ref.face_size(0) + grid_ref.face_size(1)));
	for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i <= grid_ref.Nx; i++)
	{
		int face_ind = face_x_off + grid_ref.face_ind(i, j, 0);
		Vector2i face = Vector2i(i, j);
		y[face_ind] = x(0, face);
	}

	for (int j = 0; j <= grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
	{
		int face_ind = face_y_off + grid_ref.face_ind(i, j, 1);
		Vector2i face = Vector2i(i, j);
		y[face_ind] = x(1, face);
	}
}

template<typename T>
void simplexSupport::solver2simplex(const grid2D& grid_ref, T *x, Field<T, 2>& y)
{
	const int cell_off = 0;

	for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
	{
		int cell_ind = cell_off + grid_ref.cell_ind(i, j);
		Vector2i cell = Vector2i(i, j);
		y(cell) = x[cell_ind];
	}
}

template<typename T>
void simplexSupport::solver2simplex(const grid2D& grid_ref, T *x, FaceField<T, 2>& y)
{
	const int face_x_off = 0;
	const int face_y_off = face_x_off + grid_ref.face_size(0);

	for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i <= grid_ref.Nx; i++)
	{
		int face_ind = face_x_off + grid_ref.face_ind(i, j, 0);
		Vector2i face = Vector2i(i, j);
		y(0, face) = x[face_ind];
	}

	for (int j = 0; j <= grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
	{
		int face_ind = face_y_off + grid_ref.face_ind(i, j, 1);
		Vector2i face = Vector2i(i, j);
		y(1, face) = x[face_ind];
	}
}

template<typename A, typename B, typename F>
void simplexSupport::simplex2solver(const Field<A, 2>& x, const grid2D & grid_ref, B *y, F f)
{
	const int cell_off = 0;

	memset(y + cell_off, 0, sizeof(B)*grid_ref.cell_size());

	for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
	{
		int cell_ind = cell_off + grid_ref.cell_ind(i, j);
		Vector2i cell = Vector2i(i, j);
		y[cell_ind] = f(x(cell));
	}
}

template<typename A, typename B, typename F>
void simplexSupport::simplex2solver(const FaceField<A, 2>& x, const grid2D & grid_ref, B *y, F f)
{
	const int face_x_off = 0;
	const int face_y_off = face_x_off + grid_ref.face_size(0);

	memset(y + face_x_off, 0, sizeof(B)*(grid_ref.face_size(0) + grid_ref.face_size(1)));

	for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i <= grid_ref.Nx; i++)
	{
		int face_ind = face_x_off + grid_ref.face_ind(i, j, 0);
		Vector2i face = Vector2i(i, j);
		y[face_ind] = f(x(0, face));
	}

	for (int j = 0; j <= grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
	{
		int face_ind = face_y_off + grid_ref.face_ind(i, j, 1);
		Vector2i face = Vector2i(i, j);
		y[face_ind] = f(x(1, face));
	}
}

template<typename A, typename B, typename F>
void simplexSupport::solver2simplex(const grid2D& grid_ref, A *x, Field<B, 2>& y, F f)
{
	const int cell_off = 0;

	for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
	{
		int cell_ind = cell_off + grid_ref.cell_ind(i, j);
		Vector2i cell = Vector2i(i, j);
		y(cell) = f(x[cell_ind]);
	}
}

template<typename A, typename B, typename F>
void simplexSupport::solver2simplex(const grid2D& grid_ref, A *x, FaceField<B, 2>& y, F f)
{
	const int face_x_off = 0;
	const int face_y_off = face_x_off + grid_ref.face_size(0);

	for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i <= grid_ref.Nx; i++)
	{
		int face_ind = face_x_off + grid_ref.face_ind(i, j, 0);
		Vector2i face = Vector2i(i, j);
		y(0, face) = f(x[face_ind]);
	}

	for (int j = 0; j <= grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
	{
		int face_ind = face_y_off + grid_ref.face_ind(i, j, 1);
		Vector2i face = Vector2i(i, j);
		y(1, face) = f(x[face_ind]);
	}
}

template<typename T>
void simplexSupport::simplex2solver(const Field<T, 3>& x, const grid3D & grid_ref, T * y, FieldType type)
{
	if (type == Cell)
	{
		const int cell_off = 0;

		memset(y, 0, sizeof(T)*grid_ref.cell_size());

		for (int k = 0; k < grid_ref.Nz; k++) for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
		{
			int cell_ind = cell_off + grid_ref.cell_ind(i, j, k);
			Vector3i cell = Vector3i(i, j, k);
			y[cell_ind] = x(cell);
		}
	}
	else if (type == Node)
	{
		const int node_off = 0;

		memset(y, 0, sizeof(T)*grid_ref.node_size());

		for (int k = 0; k <= grid_ref.Nz; k++) for (int j = 0; j <= grid_ref.Ny; j++) for (int i = 0; i <= grid_ref.Nx; i++)
		{
			int node_ind = node_off + grid_ref.node_ind(i, j, k);
			Vector3i node = Vector3i(i, j, k);
			y[node_ind] = x(node);
		}
	}
}

template<typename T>
void simplexSupport::simplex2solver(const EdgeField<T, 3>& x, const grid3D & grid_ref, T * y)
{
	memset(y, 0, sizeof(T)*(grid_ref.edge_size(0)+ grid_ref.edge_size(1)+ grid_ref.edge_size(2)));

	int off = 0;
	for (int d = 0; d < 3; d++)
	{
		for (int k = 0; k < grid_ref.Nz + (d != 2); k++) for (int j = 0; j < grid_ref.Ny + (d != 1); j++) for (int i = 0; i < grid_ref.Nx + (d != 0); i++)
		{
			int edge_ind = off + grid_ref.edge_ind(i, j, k, d);
			Vector3i edge = Vector3i(i, j, k);
			y[edge_ind] = x(d, edge);
		}
		off += grid_ref.edge_size(d);
	}
}

template<typename T>
void simplexSupport::simplex2solver(const FaceField<T, 3>& x, const grid3D & grid_ref, T * y)
{
	const int face_x_off = 0;
	const int face_y_off = face_x_off + grid_ref.face_size(0);
	const int face_z_off = face_y_off + grid_ref.face_size(1);

	memset(y, 0, sizeof(T)*(grid_ref.face_size(0) + grid_ref.face_size(1) + grid_ref.face_size(2)));

	for (int k = 0; k < grid_ref.Nz; k++) for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i <= grid_ref.Nx; i++)
	{
		int face_ind = face_x_off + grid_ref.face_ind(i, j, k, 0);
		Vector3i face = Vector3i(i, j, k);
		y[face_ind] = x(0, face);
	}

	for (int k = 0; k < grid_ref.Nz; k++) for (int j = 0; j <= grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
	{
		int face_ind = face_y_off + grid_ref.face_ind(i, j, k, 1);
		Vector3i face = Vector3i(i, j, k);
		y[face_ind] = x(1, face);
	}

	for (int k = 0; k <= grid_ref.Nz; k++) for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
	{
		int face_ind = face_z_off + grid_ref.face_ind(i, j, k, 2);
		Vector3i face = Vector3i(i, j, k);
		y[face_ind] = x(2, face);
	}
}

template<typename T>
void simplexSupport::solver2simplex(const grid3D & grid_ref, T * x, Field<T, 3>& y, FieldType type)
{
	if (type == Cell)
	{
		const int cell_off = 0;

		for (int k = 0; k < grid_ref.Nz; k++) for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
		{
			int cell_ind = cell_off + grid_ref.cell_ind(i, j, k);
			Vector3i cell = Vector3i(i, j, k);
			y(cell) = x[cell_ind];
		}
	}
	else if (type == Node)
	{
		const int node_off = 0;

		for (int k = 0; k <= grid_ref.Nz; k++) for (int j = 0; j <= grid_ref.Ny; j++) for (int i = 0; i <= grid_ref.Nx; i++)
		{
			int node_ind = node_off + grid_ref.node_ind(i, j, k);
			Vector3i node = Vector3i(i, j, k);
			y(node) = x[node_ind];
		}
	}
}

template<typename T>
void simplexSupport::solver2simplex(const grid3D & grid_ref, T * x, EdgeField<T, 3>& y)
{
	int off = 0;
	for (int d = 0; d < 3; d++)
	{
		for (int k = 0; k < grid_ref.Nz + (d != 2); k++) for (int j = 0; j < grid_ref.Ny + (d != 1); j++) for (int i = 0; i < grid_ref.Nx + (d != 0); i++)
		{
			int edge_ind = off + grid_ref.edge_ind(i, j, k, d);
			Vector3i edge = Vector3i(i, j, k);
			y(d, edge) = x[edge_ind];
		}
		off += grid_ref.edge_size(d);
	}
}

template<typename T>
void simplexSupport::solver2simplex(const grid3D & grid_ref, T * x, FaceField<T, 3>& y)
{
	const int face_x_off = 0;
	const int face_y_off = face_x_off + grid_ref.face_size(0);
	const int face_z_off = face_y_off + grid_ref.face_size(1);

	for (int k = 0; k < grid_ref.Nz; k++) for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i <= grid_ref.Nx; i++)
	{
		int face_ind = face_x_off + grid_ref.face_ind(i, j, k, 0);
		Vector3i face = Vector3i(i, j, k);
		y(0, face) = x[face_ind];
	}

	for (int k = 0; k < grid_ref.Nz; k++) for (int j = 0; j <= grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
	{
		int face_ind = face_y_off + grid_ref.face_ind(i, j, k, 1);
		Vector3i face = Vector3i(i, j, k);
		y(1, face) = x[face_ind];
	}

	for (int k = 0; k <= grid_ref.Nz; k++) for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
	{
		int face_ind = face_z_off + grid_ref.face_ind(i, j, k, 2);
		Vector3i face = Vector3i(i, j, k);
		y(2, face) = x[face_ind];
	}
}

template<typename A, typename B, typename F>
void simplexSupport::simplex2solver(const Field<A, 3>& x, const grid3D & grid_ref, B * y, F f, FieldType type)
{
	if (type == Cell)
	{
		const int cell_off = 0;

		memset(y, 0, sizeof(B)*grid_ref.cell_size());

		for (int k = 0; k < grid_ref.Nz; k++) for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
		{
			int cell_ind = cell_off + grid_ref.cell_ind(i, j, k);
			Vector3i cell = Vector3i(i, j, k);
			y[cell_ind] = f(x(cell));
		}
	}
	else if (type == Node)
	{
		const int node_off = 0;

		memset(y, 0, sizeof(B)*grid_ref.node_size());

		for (int k = 0; k <= grid_ref.Nz; k++) for (int j = 0; j <= grid_ref.Ny; j++) for (int i = 0; i <= grid_ref.Nx; i++)
		{
			int node_ind = node_off + grid_ref.node_ind(i, j, k);
			Vector3i node = Vector3i(i, j, k);
			y[node_ind] = f(x(node));
		}
	}
}

template<typename A, typename B, typename F>
void simplexSupport::simplex2solver(const EdgeField<A, 3>& x, const grid3D & grid_ref, B * y, F f)
{
	memset(y, 0, sizeof(B)*(grid_ref.edge_size(0) + grid_ref.edge_size(1) + grid_ref.edge_size(2)));

	int off = 0;
	for (int d = 0; d < 3; d++)
	{
		for (int k = 0; k < grid_ref.Nz + (d != 2); k++) for (int j = 0; j < grid_ref.Ny + (d != 1); j++) for (int i = 0; i < grid_ref.Nx + (d != 0); i++)
		{
			int edge_ind = off + grid_ref.edge_ind(i, j, k, d);
			Vector3i edge = Vector3i(i, j, k);
			y[edge_ind] = f(x(d, edge));
		}
		off += grid_ref.edge_size(d);
	}
}

template<typename A, typename B, typename F>
void simplexSupport::simplex2solver(const FaceField<A, 3>& x, const grid3D & grid_ref, B * y, F f)
{
	const int face_x_off = 0;
	const int face_y_off = face_x_off + grid_ref.face_size(0);
	const int face_z_off = face_y_off + grid_ref.face_size(1);

	memset(y, 0, sizeof(B)*(grid_ref.face_size(0) + grid_ref.face_size(1) + grid_ref.face_size(2)));

	for (int k = 0; k < grid_ref.Nz; k++) for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i <= grid_ref.Nx; i++)
	{
		int face_ind = face_x_off + grid_ref.face_ind(i, j, k, 0);
		Vector3i face = Vector3i(i, j, k);
		y[face_ind] = f(x(0, face));
	}

	for (int k = 0; k < grid_ref.Nz; k++) for (int j = 0; j <= grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
	{
		int face_ind = face_y_off + grid_ref.face_ind(i, j, k, 1);
		Vector3i face = Vector3i(i, j, k);
		y[face_ind] = f(x(1, face));
	}

	for (int k = 0; k <= grid_ref.Nz; k++) for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
	{
		int face_ind = face_z_off + grid_ref.face_ind(i, j, k, 2);
		Vector3i face = Vector3i(i, j, k);
		y[face_ind] = f(x(2, face));
	}

}

template<typename A, typename B, typename F>
void simplexSupport::solver2simplex(const grid3D & grid_ref, A * x, Field<B, 3>& y, F f, FieldType type)
{
	if (type == Cell)
	{
		const int cell_off = 0;

		for (int k = 0; k < grid_ref.Nz; k++) for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
		{
			int cell_ind = cell_off + grid_ref.cell_ind(i, j, k);
			Vector3i cell = Vector3i(i, j, k);
			y(cell) = f(x[cell_ind]);
		}
	}
	else if (type == Node)
	{
		const int node_off = 0;

		for (int k = 0; k <= grid_ref.Nz; k++) for (int j = 0; j <= grid_ref.Ny; j++) for (int i = 0; i <= grid_ref.Nx; i++)
		{
			int node_ind = node_off + grid_ref.node_ind(i, j, k);
			Vector3i node = Vector3i(i, j, k);
			y(node) = f(x[cell_ind]);
		}
	}
}

template<typename A, typename B, typename F>
void simplexSupport::solver2simplex(const grid3D & grid_ref, A * x, EdgeField<B, 3>& y, F f)
{
	int off = 0;
	for (int d = 0; d < 3; d++)
	{
		for (int k = 0; k < grid_ref.Nz + (d != 2); k++) for (int j = 0; j < grid_ref.Ny + (d != 1); j++) for (int i = 0; i < grid_ref.Nx + (d != 0); i++)
		{
			int edge_ind = off + grid_ref.edge_ind(i, j, k, d);
			Vector3i edge = Vector3i(i, j, k);
			y(d, edge) = f(x[edge_ind]);
		}
		off += grid_ref.edge_size(d);
	}
}

template<typename A, typename B, typename F>
void simplexSupport::solver2simplex(const grid3D & grid_ref, A * x, FaceField<B, 3>& y, F f)
{
	const int face_x_off = 0;
	const int face_y_off = face_x_off + grid_ref.face_size(0);
	const int face_z_off = face_y_off + grid_ref.face_size(1);

	for (int k = 0; k < grid_ref.Nz; k++) for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i <= grid_ref.Nx; i++)
	{
		int face_ind = face_x_off + grid_ref.face_ind(i, j, k, 0);
		Vector3i face = Vector3i(i, j, k);
		y(0, face) = f(x[face_ind]);
	}

	for (int k = 0; k < grid_ref.Nz; k++) for (int j = 0; j <= grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
	{
		int face_ind = face_y_off + grid_ref.face_ind(i, j, k, 1);
		Vector3i face = Vector3i(i, j, k);
		y(1, face) = f(x[face_ind]);
	}

	for (int k = 0; k <= grid_ref.Nz; k++) for (int j = 0; j < grid_ref.Ny; j++) for (int i = 0; i < grid_ref.Nx; i++)
	{
		int face_ind = face_z_off + grid_ref.face_ind(i, j, k, 2);
		Vector3i face = Vector3i(i, j, k);
		y(2, face) = f(x[face_ind]);
	}
}
