#include "FluidTopoOpt3D.h"
#include "Timer.h"

real FluidTopoOpt3D::Compute_Objective(const real * var)
{
	MMA2Fluid(var, design);
	return (real)Obj();
}

void FluidTopoOpt3D::Compute_Gradient(const real * var, real * grad)
{
	MMA2Fluid(var, design);
	Update_Grad();
	Fluid2MMA(design_grad, grad);

	//Numerical_Derivative();
}

void FluidTopoOpt3D::Compute_Constraint(const real * var, real * constraint)
{
	constraint[0] = 0.0;
	for (int i = 0; i < n_var; i++) {
		constraint[0] += var[i];
		//std::cout << "var[" << i << "]" << var[i];
	}
	constraint[0] = constraint[0] / ((real)n_var) - frac;
	constraint[0] *= constraint_coef;
	//constraint[0]=pow((constraint[0]/((real)n_var)-frac),2);
	std::cout << "constraint: " << constraint[0] << std::endl;
}

void FluidTopoOpt3D::Compute_Constraint_Grad(const real * var, real * constraint_grad)
{
	for (int i = 0; i < n_var; i++) {
		constraint_grad[i * n_cons] = constraint_coef * 1.0 / (real)n_var;
		//constraint_grad[i*n_cons]=2.0/((real)n_var*(real)n_var)*var[i]-2.0*frac/((real)n_var);
	}
	//std::cout << "************" << constraint_grad[0] << std::endl;
}

void FluidTopoOpt3D::Write_Substep(const int frame)
{
	std::cout << "Optimization iteration " << frame << std::endl;
	if (Write_Output_Files_Callback)Write_Output_Files_Callback(frame);
}

void FluidTopoOpt3D::Init(int _grid_size, Scalar _frac, int _max_iter)
{
	grid_size = _grid_size;

	// no dx in solver, lower alpha_os instead
	alpha_os /= pow((Scalar)grid_size / 32, 1);

	frac = _frac;
	max_iter_num = _max_iter;

	Vector3i sz = Vector3i(grid_size, grid_size, grid_size);

	grid.Initialize(sz, 1.);
	mac_grid.Initialize(grid);

	steady.init(grid_size);
	energy.init(grid_size);

	cell_penalty.Resize(sz);
	face_penalty.Resize(sz);

	design.Resize(sz);
	design_grad.Resize(sz);

	buffer.Resize(sz);

	Init_MMA();
}

void FluidTopoOpt3D::Init_MMA()
{
	n_var = grid.Number_Of_Cells();
	n_cons = 1;
	//max_iter_num = 1000;
	movlim = (real).2;
	var_lb = (real).01;
	var_ub = (real)1.;
	Allocate_Data();
}

void FluidTopoOpt3D::Init_Boundary(const Field<Scalar, 3>& cell_vol, const Field<int, 3>& cell_fixed, const Field<Scalar, 3>& _cell_b, \
	const FaceField<Scalar, 3>& face_vol, const FaceField<int, 3>& face_fixed, const FaceField<Scalar, 3>& _face_b, const Field<Scalar, 3>& _cell_penalty)
{
	steady.init_boundary(cell_vol, cell_fixed, face_vol, face_fixed, closed);
	energy.init_boundary(cell_vol, cell_fixed, face_vol, face_fixed);

	cell_b = _cell_b; face_b = _face_b; cell_penalty = _cell_penalty;

	//iterate_cell_d(iter, grid, 2)
	//{
	//	const Vector2i& cell = iter.Coord();
	//	design(cell) = cell_fixed(cell) == 0 ? frac : (Scalar)0;
	//}

	iterate_cell_d(iter, grid, 3)
	{
		const Vector3i& cell = iter.Coord();
		design(cell) = frac;
	}

	Fluid2MMA(design, var);
}

Scalar FluidTopoOpt3D::alpha(Scalar rho)
{
	//return rho;
	return alpha_os + (alpha_us - alpha_os) * rho * ((Scalar)1 + q) / (rho + q);
	//return -(alpha_os + (alpha_us - alpha_os) * rho * ((Scalar)1 + q) / (rho + q));
}

Scalar FluidTopoOpt3D::D_alpha(Scalar rho)
{
	//return (Scalar)1;
	return (alpha_us - alpha_os) * (((Scalar)1 + q) / (rho + q) - rho * ((Scalar)1 + q) / pow(rho + q, 2));
	//return -((alpha_us - alpha_os) * (((Scalar)1 + q) / (rho + q) - rho * ((Scalar)1 + q) / pow(rho + q, 2)));

}

void FluidTopoOpt3D::MMA2Fluid(const real *v, Field<Scalar, 3>& f)
{
	std::memcpy(&buffer.array[0], v, sizeof(real)*n_var);
	iterate_cell_d(iter, grid, 3)
	{
		const Vector3i& cell = iter.Coord();
		f(cell) = (Scalar)buffer(cell);
	}
}

void FluidTopoOpt3D::Fluid2MMA(const Field<Scalar, 3> f, real *v)
{
	iterate_cell_d(iter, grid, 3)
	{
		const Vector3i& cell = iter.Coord();
		buffer(cell) = (real)f(cell);
	}
	std::memcpy(v, &buffer.array[0], sizeof(real)*n_var);
}

void FluidTopoOpt3D::Update_Fluid()
{
	Timer<real> timer;
	timer.Reset();
	face_penalty.Fill((Scalar)0);
	iterate_cell_d(iter, grid, 3)
	{
		const Vector3i& cell = iter.Coord();
		Scalar t = (Scalar).5*alpha(pow(design(cell), 3));
		for (int dd = 0; dd < 3; dd++)
		{
			face_penalty(dd, cell) += t;
			face_penalty(dd, cell + Vector3i::Unit(dd)) += t;
		}
	}
	steady.update_penalty(cell_penalty, face_penalty);
	steady.update_b(cell_b, face_b);
	steady.solve();

	vel = steady.vel;
	p = steady.p;
	energy.update_penalty(cell_penalty, face_penalty);
	timer.Elapse_And_Output("steady state");	
}

void FluidTopoOpt3D::Update_Grad()
{
	Timer<real> timer;
	timer.Reset();
	//Update_Fluid();			// called in obj
	energy.compute_gradient_x(vel, p);
	energy.compute_gradient_q(vel, p);
	steady.update_b(energy.p_grad, energy.vel_grad);
	steady.solve();

	temp_grad = vel;
	temp_grad *= steady.vel;
	temp_grad *= (Scalar)-1;
	temp_grad += energy.face_penalty_grad;

	iterate_cell_d(iter, grid, 3)
	{
		const Vector3i& cell = iter.Coord();
		Scalar t = (Scalar).5*D_alpha(pow(design(cell), 3))*(Scalar)3 * pow(design(cell), 2);
		Scalar sum = (Scalar)0;
		for (int dd = 0; dd < 3; dd++)
		{
			sum += temp_grad(dd, cell);
			sum += temp_grad(dd, cell+Vector3i::Unit(dd));
		}
		design_grad(cell) = sum * t;
	}
	timer.Elapse_And_Output("gradient");
}

Scalar FluidTopoOpt3D::Obj()
{
	Update_Fluid();
	return energy.compute_energy(vel, p);
}

void FluidTopoOpt3D::Numerical_Derivative()
{
	Field nd = design;
	Scalar dx = 1e-2;
	for (int i = 0; i < design.array.size(); i++)
	{
		Scalar origin = design(i);
		design(i) = origin + dx;
		Scalar obj1 = Obj();
		design(i) = origin - dx;
		Scalar obj2 = Obj();
		nd(i) = (obj1 - obj2) / ((Scalar)2 * dx);
		design(i) = origin;
	}
	Update_Grad();
	for (int i = 0; i < design.array.size(); i++)
		printf("cell:%d ana:%lf num:%lf\n", i, design_grad(i), nd(i));
}
