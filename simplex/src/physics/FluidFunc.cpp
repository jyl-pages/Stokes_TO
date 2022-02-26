#include "FluidFunc.h"
#include "BoundaryCondition.h"

namespace FluidFunc {

	//////////////////////////////////////////////////////////////////////////
	////vorticity helper functions
	//////store vorticity on cell center
	//template<int dim> real dudx_on_cell(int u_axis,int x_axis,const Vector<int,dim>& cell,const MacGrid<dim>& mac_grid,const Interpolation<dim>& intp,const FaceField<real,dim>& v);
	//////store vorticity on node
	//template<int dim> real dudx_on_node(int u_axis,int x_axis,const Vector<int,dim>& node,const MacGrid<dim>& mac_grid,const FaceField<real,dim>& v);

	template<int dim> real dudx_on_cell(int u_axis,int x_axis,const Vector<int,dim>& cell,const MacGrid<dim>& mac_grid,const Interpolation<dim>& intp,const FaceField<real,dim>& v)
	{

		Vector<int,dim> nb1=cell-Vector<int,dim>::Unit(x_axis);
		real u1=(real)0;if(mac_grid.grid.Valid_Cell(nb1))u1=intp.Interpolate_Faces_To_Cell(v,nb1,u_axis);
		Vector<int,dim> nb2=cell+Vector<int,dim>::Unit(x_axis);
		real u2=(real)0;if(mac_grid.grid.Valid_Cell(nb2))u2=intp.Interpolate_Faces_To_Cell(v,nb2,u_axis);
		return u2-u1;
	}

	template<int dim> real dudx_on_node(int u_axis,int x_axis,const Vector<int,dim>& node,const MacGrid<dim>& mac_grid,const FaceField<real,dim>& v)
	{
		Vector<int,dim> nb1=node-Vector<int,dim>::Unit(x_axis);
		real u1=(real)0;if(mac_grid.Valid_Face(u_axis,nb1))u1=v(u_axis,nb1);
		Vector<int,dim> nb2=node;
		real u2=(real)0;if(mac_grid.Valid_Face(u_axis,nb2))u2=v(u_axis,nb2);
		return u2-u1;
	}

	void Curl_On_Cell(const MacGrid<2>& mac_grid,const FaceField<real,2>& v,Field<real,2>& vorticity)
	{
		Interpolation<2> intp(mac_grid);
		int v_axis[2]={ 1,0 }; real coef_axis[2]={ 1.,-1. };
		real one_over_two_dx=(real).5 / mac_grid.grid.dx;
		//iterate_cell(iter,mac_grid.grid){const VectorDi& cell=iter.Coord();
		int cell_num=mac_grid.grid.cell_counts.prod();
		#pragma omp parallel for
		for (int i=0; i<cell_num; i++) {
			Vector<int,2> cell=mac_grid.grid.Cell_Coord(i);
			vorticity(cell)=(dudx_on_cell<2>(1,0,cell,mac_grid,intp,v)-dudx_on_cell<2>(0,1,cell,mac_grid,intp,v)) * one_over_two_dx;
		}	////dvdx-dudy
	}

	void Curl_On_Cell(const MacGrid<3>& mac_grid,const FaceField<real,3>& v,Field<Vector3,3>& vorticity) {
		Interpolation<3> intp(mac_grid);
		real one_over_two_dx=(real).5 / mac_grid.grid.dx;
		int cell_num=mac_grid.grid.cell_counts.prod();
		#pragma omp parallel for
		for (int i=0; i<cell_num; i++) {
			Vector<int,3> cell=mac_grid.grid.Cell_Coord(i);
			Vector3 vor;	////(dwdy-dvdz,dudz-dwdx,dvdx-dudy)
			vor[0]=(dudx_on_cell<3>(2,1,cell,mac_grid,intp,v)-dudx_on_cell<3>(1,2,cell,mac_grid,intp,v)) * one_over_two_dx;
			vor[1]=(dudx_on_cell<3>(0,2,cell,mac_grid,intp,v)-dudx_on_cell<3>(2,0,cell,mac_grid,intp,v)) * one_over_two_dx;
			vor[2]=(dudx_on_cell<3>(1,0,cell,mac_grid,intp,v)-dudx_on_cell<3>(0,1,cell,mac_grid,intp,v)) * one_over_two_dx;
			vorticity(cell)=vor;}
	}

	void Curl_On_Node(const MacGrid<2>& mac_grid,const FaceField<real,2>& v,Field<real,2>& vorticity)
	{
		real one_over_dx=(real)1 / mac_grid.grid.dx;
		int node_num=mac_grid.grid.Number_Of_Nodes();
		#pragma omp parallel for
		for (int i=0; i<node_num; i++) {
			Vector<int,2> node=mac_grid.grid.Node_Coord(i);
			vorticity(node)=(dudx_on_node<2>(1,0,node,mac_grid,v)-dudx_on_node<2>(0,1,node,mac_grid,v)) * one_over_dx;
		}	////dvdx-dudy
	}

	void Curl_On_Face(const MacGrid<2>& mac_grid,Field<real,2>& psi,FaceField<real,2>& vorticity)
	{
		const int d=2; Typedef_VectorDii(d);
		real one_over_dx=(real)1 / mac_grid.grid.dx;
		iterate_face_in_one_dim(0,iter,mac_grid) {
			const VectorDi& face=iter.Coord();
			vorticity(0,face)=(psi(face+VectorDi::Unit(1))-psi(face)) * one_over_dx;}
		iterate_face_in_one_dim(1,iter,mac_grid) {
			const VectorDi& face=iter.Coord();
			vorticity(1,face)=-(psi(face+VectorDi::Unit(0))-psi(face)) * one_over_dx;}
	}

	template<int d>
	void Enforce_Vorticity_Confinement(const MacGrid<d>& mac_grid, const real dt, FaceField<real, d>& velocity, const real vor_conf_coef)
	{
		Typedef_VectorDii(d);
		std::conditional_t<d == 2, Field<real, 2>, Field<VectorD, 3> > vorticity;
		vorticity.Resize(mac_grid.grid.cell_counts);
		Curl_On_Cell(mac_grid, velocity, vorticity);
		Field<VectorD, d> vorticity_confinement_force(mac_grid.grid.cell_counts, VectorD::Zero());

		//N = (grad(|vor|)) / |grad(|vor|)|
		int cell_num = mac_grid.grid.Number_Of_Cells();
		real dx = mac_grid.grid.dx;
		real one_over_2dx = (real)0.5 / dx;
#pragma omp parallel for
		for (int i = 0; i < cell_num; i++) {
			const VectorDi& cell = mac_grid.grid.Cell_Coord(i);
			if (mac_grid.grid.Is_Boundary_Cell(cell)) { continue; }	//Can't calculate vorticity norm gradient if cell is on boundary
			VectorD pos = mac_grid.grid.Center(cell);
			VectorD vor_norm_grad; //vorticity norm gradient direction
			for (int i = 0; i < d; i++) {
				VectorDi cell_left = cell - VectorDi::Unit(i); VectorDi cell_right = cell + VectorDi::Unit(i);
				real vorticity_left_norm; real vorticity_right_norm;
				if constexpr (d == 2) { vorticity_left_norm = abs(vorticity(cell_left)); vorticity_right_norm = abs(vorticity(cell_right)); }
				else { vorticity_left_norm = vorticity(cell_left).norm(); vorticity_right_norm = vorticity(cell_right).norm(); }
				vor_norm_grad[i] = (vorticity_right_norm - vorticity_left_norm) * one_over_2dx;
			}
			vor_norm_grad.normalize();
			if constexpr (d == 2) { vorticity_confinement_force(cell) = vor_conf_coef * dx * (-AuxFunc::Cross(vorticity(cell), vor_norm_grad)); }
			else { vorticity_confinement_force(cell) = vor_conf_coef * dx * vor_norm_grad.cross(vorticity(cell)); }
		}

		//Apply vorticity confinement force to velocity
		Interpolation<d> intp(mac_grid.grid);
		for (int axis = 0; axis < d; axis++) {
			int face_num = mac_grid.Number_Of_Faces(axis);
#pragma omp parallel for
			for (int i = 0; i < face_num; i++) {
				VectorDi face = mac_grid.Face_Coord(axis, i);
				VectorD pos = mac_grid.Face_Center(axis, face);
				velocity(axis, face) += dt * intp.Interpolate_Centers(vorticity_confinement_force, pos)[axis];
			}
		}
	}
	template void Enforce_Vorticity_Confinement<2>(const MacGrid<2>& mac_grid, const real dt, FaceField<real, 2>& velocity, const real vor_conf_coef);
	template void Enforce_Vorticity_Confinement<3>(const MacGrid<3>& mac_grid, const real dt, FaceField<real, 3>& velocity, const real vor_conf_coef);

	template<int d>
	void Correct_Velocity_With_Pressure(const MacGrid<d>& mac_grid, FaceField<real, d>& velocity, const Field<real, d>& pressure, const FaceField<real, d>& alpha, const Field<ushort, d>& type)
	{
		Typedef_VectorDii(d);
		for (int axis = 0; axis < d; axis++) {
			int face_num = mac_grid.face_grids[axis].node_counts.prod();
#pragma omp parallel for
			for (int i = 0; i < face_num; i++) {
				VectorDi face = mac_grid.face_grids[axis].Node_Coord(i);
				VectorDi cell[2]; for (int i = 0; i < 2; i++)cell[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
				real cell_p[2]; for (int i = 0; i < 2; i++) {
					if (mac_grid.grid.Valid_Cell(cell[i]) && type(cell[i])== (ushort)CellType::Fluid) {
						cell_p[i] = pressure(cell[i]);
					}
					else cell_p[i] = 0;
				}
				real gradp = (cell_p[1] - cell_p[0]);
				velocity(axis, face) -= alpha(axis, face) * (cell_p[1] - cell_p[0]);
			}
		}
	}
	template void Correct_Velocity_With_Pressure<2>(const MacGrid<2>& mac_grid, FaceField<real, 2>& velocity, const Field<real, 2>& pressure, const FaceField<real, 2>& alpha, const Field<ushort, 2>& type);
	template void Correct_Velocity_With_Pressure<3>(const MacGrid<3>& mac_grid, FaceField<real, 3>& velocity, const Field<real, 3>& pressure, const FaceField<real, 3>& alpha, const Field<ushort, 3>& type);

	//////////////////////////////////////////////////////////////////////////
	////diffusion

	template<class T,int d> void Diffusion(Field<T,d>& field,const Grid<d>& grid,const real a,std::function<real(const Vector<real,d>&)> Phi/*=nullptr*/,const int iter_num/*=10*/)
	{Typedef_VectorDii(d);
		if(Phi==nullptr){	////diffusion without boundary
			int cell_num=field.counts.prod();
			Field<T,d> field_0=field;
			for(int k=0;k<iter_num;k++)for(int i=0;i<cell_num;i++){
				VectorDi cell=grid.Cell_Coord(i);
				T delta=Zero<T>();int nb_n=0;
				for(int j=0;j<grid.Number_Of_Nb_C();j++){
					VectorDi nb_cell=grid.Nb_C(cell,j);
					if(!grid.Valid_Cell(nb_cell))continue;
					delta+=field_0(nb_cell);nb_n++;}
				field(cell)=(field_0(cell)+a*delta)/((real)1+(nb_n)*a);}}
		else{	////diffusion with boundary
			int cell_num=field.counts.prod();
				Field<T,d> field_0=field;
				for(int k=0;k<iter_num;k++)for(int i=0;i<cell_num;i++){
					VectorDi cell=grid.Cell_Coord(i);
					VectorD xi=grid.Center(cell);
					real phi_i=Phi(xi);if(phi_i>=(real)0)continue;

					T delta=Zero<T>();int nb_n=0;
					for(int j=0;j<grid.Number_Of_Nb_C();j++){
						VectorDi nb_cell=grid.Nb_C(cell,j);
						if(!grid.Valid_Cell(nb_cell))continue;
						VectorD xj=grid.Center(nb_cell);
						real phi_j=Phi(xj);if(phi_j>=(real)0)continue;
						delta+=field_0(nb_cell);nb_n++;}
					field(cell)=(field_0(cell)+a*delta)/((real)1+(nb_n)*a);}}
	}

	template void Diffusion<real,2>(Field<real,2>&,const Grid<2>&,const real,std::function<real(const Vector<real,2>&)>,const int);
	template void Diffusion<real,3>(Field<real,3>&,const Grid<3>&,const real,std::function<real(const Vector<real,3>&)>,const int);
	template void Diffusion<Vector2,2>(Field<Vector2,2>&,const Grid<2>&,const real,std::function<real(const Vector<real,2>&)>,const int);
	template void Diffusion<Vector3,3>(Field<Vector3,3>&,const Grid<3>&,const real,std::function<real(const Vector<real,3>&)>,const int);
	//template void Diffusion<Vector4,2>(Field<Vector4,2>&,const Grid<2>&,const real,std::function<real(const Vector<real,2>&)>,const int);
	//template void Diffusion<Vector4,3>(Field<Vector4,3>&,const Grid<3>&,const real,std::function<real(const Vector<real,3>&)>,const int);

	template<int d> void Diffusion(FaceField<real,d>& field,const MacGrid<d>& mac_grid,const real a,std::function<real(const Vector<real,d>&)> Phi/*=nullptr*/,const int iter_num/*=10*/)
	{Typedef_VectorDii(d);
		////TODO: parallellization
		if(Phi==nullptr){	////diffusion without boundary
			FaceField<real,d> field_0=field;
			for(int k=0;k<iter_num;k++){
				iterate_face(axis,iter,mac_grid){
					const VectorDi& face=iter.Coord();
					real delta=(real)0;int nb_n=0;
					for(int j=0;j<mac_grid.face_grids[axis].Number_Of_Nb_C();j++){
						VectorDi nb_face=mac_grid.face_grids[axis].Nb_C(face,j);
						if(!mac_grid.Valid_Face(axis,nb_face))continue;
						delta+=field(axis,nb_face);nb_n++;}
					field(axis,face)=(field_0(axis,face)+a*delta)/((real)1+(nb_n)*a);}}}
		else{	////diffusion with boundary
			FaceField<real,d> field_0=field;
			for(int k=0;k<iter_num;k++){
				iterate_face(axis,iter,mac_grid){
					const VectorDi& face=iter.Coord();
					VectorD xi=mac_grid.Face_Center(axis,face);
					real phi_i=Phi(xi);if(phi_i>=(real)0)continue;

					real delta=(real)0;int nb_n=0;
					for(int j=0;j<mac_grid.face_grids[axis].Number_Of_Nb_C();j++){
						VectorDi nb_face=mac_grid.face_grids[axis].Nb_C(face,j);
						if(!mac_grid.Valid_Face(axis,nb_face))continue;
						VectorD xj=mac_grid.Face_Center(axis,nb_face);
						real phi_j=Phi(xj);if(phi_j>=(real)0)continue;
						delta+=field(axis,nb_face);nb_n++;}
					field(axis,face)=(field_0(axis,face)+a*delta)/((real)1+(nb_n)*a);}}}
	}

	template void Diffusion<2>(FaceField<real,2>&,const MacGrid<2>&,const real,std::function<real(const Vector<real,2>&)>,const int);
	template void Diffusion<3>(FaceField<real,3>&,const MacGrid<3>&,const real,std::function<real(const Vector<real,3>&)>,const int);

	//////////////////////////////////////////////////////////////////////////
	////Jacobian
	
	template<int d> void Jacobian(const FaceField<real,d>& field,const Interpolation<d>& intp,const Vector<real,d>& pos,Matrix<real,d>& jac,std::function<real(const Vector<real,d>&)> Phi/*=nullptr*/)
	{Typedef_VectorDii(d);Typedef_MatrixD(d);
		if(Phi==nullptr){	////without boundary
			real dx=intp.mac_grid.grid.dx;
			for(int i=0;i<d;i++){
				VectorD p0=pos-VectorD::Unit(i)*dx;
				VectorD p1=pos+VectorD::Unit(i)*dx;
				VectorD v0=intp.Interpolate_Face_Vectors(field,p0);
				VectorD v1=intp.Interpolate_Face_Vectors(field,p1);
				jac.col(i)=(v1-v0)/((real)2*dx);}}
		else{	////with boundary
			real dx=intp.mac_grid.grid.dx;
			////assuming pos is inside phi
			for(int i=0;i<d;i++){
				bool i0=true,i1=true;
				VectorD p0=pos-VectorD::Unit(i)*dx;
				if(Phi(p0)>=0){i0=false;p0=pos;}
				VectorD p1=pos+VectorD::Unit(i)*dx;
				if(Phi(p1)>=0){i1=false;p1=pos;}
				real len=(real)2*dx;
				VectorD v0=intp.Interpolate_Face_Vectors(field,p0);
				VectorD v1=intp.Interpolate_Face_Vectors(field,p1);

				if(i0&&i1){
					VectorD v0=intp.Interpolate_Face_Vectors(field,p0);
					VectorD v1=intp.Interpolate_Face_Vectors(field,p1);
					jac.col(i)=(v1-v0)/((real)2*dx);}
				else if((!i0&&i1)||(!i1&&i0)){
					VectorD v0=intp.Interpolate_Face_Vectors(field,p0);
					VectorD v1=intp.Interpolate_Face_Vectors(field,p1);
					jac.col(i)=(v1-v0)/(dx);}
				else{jac.col(i)=VectorD::Zero();}}}
	}

	template void Jacobian<2>(const FaceField<real,2>&,const Interpolation<2>&,const Vector<real,2>&,Matrix<real,2>&,std::function<real(const Vector<real,2>&)> Phi);
	template void Jacobian<3>(const FaceField<real,3>&,const Interpolation<3>&,const Vector<real,3>&,Matrix<real,3>&,std::function<real(const Vector<real,3>&)> Phi);

	//////////////////////////////////////////////////////////////////////////
	////Deformation gradient


}