#include "CPXFunc.h"
#include "TypeFunc.h"
#include "AuxFunc.h"
//#include "FluidFunc.h"

namespace CPXFunc {
//	template<int d>
//	void Eliminate_Divergence(const MacGrid<d>& mac_grid, FaceField<real, d>& velocity, const FaceField<real, d>& alpha, const Field<ushort, d>& type, const BoundaryConditionMacGrid<d>& bc)
//	{
//		AuxFunc::Crash_With_Info("CPXFunc need to remove dependence on fluidfunc ");
//		int n = mac_grid.grid.cell_counts[0];
//		for (int axis = 0; axis < d; axis++) { if (mac_grid.grid.cell_counts[axis] != n) { AuxFunc::Crash_With_Info("Eliminate_Divergence: cell size not equal"); } }
//		if (n % 8) AuxFunc::Crash_With_Info("Eliminate_Divergence: cell size must be divided by 8");
//
//		Typedef_VectorDii(d);
//		using Projection = typename If<d == 2, Poisson, Poisson3D >::Type;
//		static Projection poisson;
//		static bool inited = false;
//		if (!inited) {
//			poisson.init(mac_grid.grid.cell_counts[0], mac_grid.grid.dx);
//			poisson.cg.verbose = false;
//			inited = true;
//		}
//		//face_vol
//		FaceField<Scalar, d> face_vol = alpha;
//		//generate cell_fixed
//		int cell_num = mac_grid.grid.cell_counts.prod();
//		Field<int, d> cell_fixed(mac_grid.grid.cell_counts, 0);
//#pragma omp parallel for
//		for (int idx = 0; idx < cell_num; idx++) {
//			const VectorDi& cell = mac_grid.grid.Cell_Coord(idx);
//			cell_fixed(cell) = (type(cell) == (ushort)CellType::Fluid) ? 0 : 1;
//		}
//		//generate cell_b
//		Field<Scalar, d> cell_b(mac_grid.grid.cell_counts, 0);
//#pragma omp parallel for
//		for (int idx = 0; idx < cell_num; idx++) {
//			const VectorDi& cell = mac_grid.grid.Cell_Coord(idx);
//			real div = 0;
//			for (int axis = 0; axis < d; axis++) {
//				VectorDi faces[2]; for (int i = 0; i < 2; i++) faces[i] = cell + VectorDi::Unit(axis) * i;
//				real vels[2];
//				for (int i = 0; i < 2; i++) {
//					//vels[i] = velocity(axis, faces[i]);
//					if (bc.Is_Psi_N(axis, faces[i])) {
//						vels[i] = bc.Psi_N_Value(axis, faces[i]);
//					}
//					else {
//						vels[i] = alpha(axis, faces[i]) * velocity(axis, faces[i]);
//					}
//				}
//				div += vels[1] - vels[0];
//				//= cell + VectorDi::Unit(axis);
//				//if(psi_)
//				//div += velocity(axis, cell + VectorDi::Unit(axis))
//				//	- velocity(axis, cell);
//			}
//			cell_b(cell) = -div;
//		}
//		//solve system
//		poisson.init_boundary(face_vol, cell_fixed, false);
//		poisson.update_b(cell_b);
//		poisson.solve();
//		Field<real, d> pressure = poisson.x;
//		//FluidFunc::Correct_Velocity_With_Pressure(mac_grid, velocity, pressure, alpha, type);
//	}
//	template void Eliminate_Divergence<2>(const MacGrid<2>& mac_grid, FaceField<real, 2>& velocity, const FaceField<real, 2>& alpha, const Field<ushort, 2>& type, const BoundaryConditionMacGrid<2>& bc);
//	template void Eliminate_Divergence<3>(const MacGrid<3>& mac_grid, FaceField<real, 3>& velocity, const FaceField<real, 3>& alpha, const Field<ushort, 3>& type, const BoundaryConditionMacGrid<3>& bc);
}