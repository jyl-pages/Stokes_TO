#pragma once

#include "OptimizerMma.h"
#include "FluidSteadyState.h"
#include "FluidEnergy.h"

#include "MacGrid.h"
#include "Grid.h"
#include "Field.h"
#include "FaceField.h"

class FluidTopoOpt : public OptimizerMMA
{
public:
	int grid_size;
	Scalar frac;

	MacGrid<2> mac_grid;
	Grid<2> grid;

	FluidSteadyState steady;
	FluidEnergy energy;

	FaceField<Scalar, 2> face_penalty;
	Field<Scalar, 2> cell_penalty;

	FaceField<Scalar, 2> face_b;
	Field<Scalar, 2> cell_b;

	Field<Scalar, 2> design;
	Field<Scalar, 2> design_grad;

	FaceField<Scalar, 2> temp_grad;

	FaceField<Scalar, 2> vel;
	Field<Scalar, 2> p;

	Field<real, 2> buffer;
	
	std::function<void(const int)> Write_Output_Files_Callback;

	// tune parameters
	Scalar q = 0.014; 

	Scalar alpha_os = 5e2;
	Scalar alpha_us = 1e-3;
	Scalar constraint_coef = 1e4;

	bool closed;

public:
	virtual real Compute_Objective(const real* var) override;

	virtual void Compute_Gradient(const real* var, real* grad) override;

	virtual void Compute_Constraint(const real* var, real* constraint) override;

	virtual void Compute_Constraint_Grad(const real* var, real* constraint_grad) override;

	virtual void Write_Substep(const int frame) override;

public:
	void Init(int _grid_size, Scalar _frac, int _max_iter);

	void Init_MMA();

	void Init_Boundary(const Field<Scalar, 2>& cell_vol, const Field<int, 2>& cell_fixed, const Field<Scalar, 2>& _cell_b, \
		const FaceField<Scalar, 2>& face_vol, const FaceField<int, 2>& face_fixed, const FaceField<Scalar, 2>& _face_b, \
		const Field<Scalar,2>& _cell_penalty);

	Scalar alpha(Scalar rho);

	Scalar D_alpha(Scalar rho);

	void MMA2Fluid(const real *v, Field<Scalar, 2>& f);

	void Fluid2MMA(const Field<Scalar, 2> f, real *v);

	void Update_Fluid();

	void Update_Grad();

	Scalar Obj();

	void Numerical_Derivative();
};
