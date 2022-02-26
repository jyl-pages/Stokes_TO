//////////////////////////////////////////////////////////////////////////
// Assemble data for opengl_viewer. Especially for Points.
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "RenderFunc.h"
#include "AuxFunc.h"
#include "Interpolation.h"

namespace RenderFunc {
	template<int d>
	void Write_Face_Scalar_Field(const std::string& file_name, const MacGrid<d>& mac_grid, const FaceField<real, d>& val)
	{
		Field<Vector<real, d>, d> grid_v;
		Field<real, d> grid_a;
		grid_v.Resize(mac_grid.grid.cell_counts, Vector<real, d>::Zero());
		grid_a.Resize(mac_grid.grid.cell_counts, 0.0);

		Interpolation<d> intp(mac_grid);
		intp.Interpolate_Faces_To_Cells(val, grid_v);
		iterate_cell(iter, mac_grid.grid) {
			Vector<int, d> cell = iter.Coord();
			Vector<real, d> cell_val = grid_v(cell);
			real a = 0; for (int k = 0; k < d; k++) a += cell_val[k];
			grid_a(cell) = a / d;
		}
		grid_a.Write_To_File_3d(file_name);
	}
	template void Write_Face_Scalar_Field<2>(const std::string& file_name, const MacGrid<2>& mac_grid, const FaceField<real, 2>& val);
	template void Write_Face_Scalar_Field<3>(const std::string& file_name, const MacGrid<3>& mac_grid, const FaceField<real, 3>& val);

	template<int d>
	void Write_Customized_Segments(std::string file_name, const Array<Vector<real, d>>& xs, const Array<Vector<real, d>>& normals, const Array<real>& len_arr, const real scale)
	{
		Array<Vector<real, d> > to_write;
		int pn = (int)len_arr.size();
		to_write.resize(pn);
#pragma omp parallel for
		for (int i = 0; i < pn; i++) {
			to_write[i] = normals[i].normalized() * len_arr[i] * scale;
		}
		Write_Segments_To_File_3d_Fast<d, real>(xs, to_write, file_name);//in Particles.h
	}
	template void Write_Customized_Segments<2>(std::string file_name, const Array<Vector2>& xs, const Array<Vector2>& normals, const Array<real>& arr, const real scale);
	template void Write_Customized_Segments<3>(std::string file_name, const Array<Vector3>& xs, const Array<Vector3>& normals, const Array<real>& arr, const real scale);

	template<int d>
	void Write_Scalars_As_Points(std::string file_name, const GeometryParticles<d>& points, const Array<real>& arr, const real scale)
	{
		int pn = points.Size();
		if (arr.size() != pn) AuxFunc::Crash_With_Info("RenderFunc::Write_Scalars_As_Points error: size not match");
		Array<Vector<real, d>> to_write(pn);
#pragma omp parallel for
		for (int i = 0; i < pn; i++) {
			to_write[i] = points.X(i) + points.Normal(i).normalized() * arr[i] * scale;
		}
		Write_To_File_3d_Fast<d, real>(to_write, file_name);
	}
	template void Write_Scalars_As_Points<2>(std::string file_name, const GeometryParticles<2>& points, const Array<real>& arr, const real scale);
	template void Write_Scalars_As_Points<3>(std::string file_name, const GeometryParticles<3>& points, const Array<real>& arr, const real scale);
}
