#include "Poisson3D.h"

#include "PoissonMapping3D.h"
//#include "PoissonDiagMapping.h"
//#include "PoissonSmoothing.h"
#include "PoissonLike.h"
#include "PoissonSubSystem3D.h"
//#include "PoissonDirectSolve.h"
#include "MultigridSolve.h"
#include "DescentSmoother.h"

#include "SimplexSupport.h"

void Poisson3D::init(Vector3i _grid_size, Scalar _dx)
{
	grid_size = _grid_size;
	grid.Initialize(grid_size, _dx);
	int n = grid_size.minCoeff();
	l = round(std::log2(n & -n)) - 1;

	x.Resize(grid_size);
	b.Resize(grid_size);

	PoissonDescriptor<3> descr = PoissonDescriptor<3>();
	descr.init(grid_size);

	temp_x = new Scalar[descr.size];
	temp_b = new Scalar[descr.size];

	cg = ConjugatedGradient();
	cg.CG_relative_residual_thres = 1e-3;
	//cg.linear_mapping = new PoissonDiagMapping(mapping);

	MultigridSolve *mg_solve = new MultigridSolve();
	mg_solve->init(l);

	mg_descr = new PoissonDescriptor<3>[l];
	mg_descr[l - 1] = descr;

	for (int i = l - 2; i >= 0; i--) mg_descr[i] = createSubSystem(mg_descr[i + 1]);
	for (int i = 0; i < l; i++)
	{
		PoissonMapping3DFixed *t0 = new PoissonMapping3DFixed();
		t0->init(&mg_descr[i]);
		mg_solve->mapping[i] = t0;
		if (i != 0)
		{
			//PoissonSmoothing *t1 = new PoissonSmoothing(t0);
			LinearMapping *t1 = PoissonLike::GenerateJacobiPreconditioner(mg_descr[i].grid, t0);
			mg_solve->preSmoother[i] = t1;
			//mg_solve->preSmoother[i] = nullptr;

			//PoissonSmoothing *t2 = new PoissonSmoothing(t0);
			mg_solve->postSmoother[i] = t1;
			//mg_solve->postSmoother[i] = nullptr;

			PoissonDownSample3D *t3 = new PoissonDownSample3D(&mg_descr[i], &mg_descr[i - 1]);
			mg_solve->downSample[i - 1] = t3;

			PoissonUpSample3D *t4 = new PoissonUpSample3D(&mg_descr[i], &mg_descr[i - 1]);
			mg_solve->upSample[i - 1] = t4;
		}
		else
		{
			//PoissonDirectSolve *direct = new PoissonDirectSolve(t0);
			//direct->finish();
			LinearMapping *direct = PoissonLike::GenerateDirectSolve(mg_descr[i].grid, t0, closed);
			mg_solve->preSmoother[i] = direct;
			mg_solve->postSmoother[i] = nullptr;
		}
	}

	mg_solve->finish();

	//cg.preconditioner = new DescentSmoother(PoissonLike::GenerateJacobiPreconditioner(descr.grid, cg.linear_mapping), 20, cg.linear_mapping);
	cg.preconditioner = mg_solve;
	cg.linear_mapping = mg_solve->mapping[l - 1];
	cg.Init(50);
}

void Poisson3D::solve()
{
	//PoissonSmoothing *temp = new PoissonSmoothing((PoissonMapping*)cg.linear_mapping);
	for (int i = l - 2; i >= 0; i--) updateSubSystem(mg_descr[i], mg_descr[i + 1]);
	//if (cg.preconditioner) cg.preconditioner = nullptr;
	cg.preconditioner->update();

	cudaDeviceSynchronize();
	checkCudaErrors(cudaGetLastError());

	const PoissonDescriptor<3>& descr = mg_descr[l - 1];

	simplexSupport::simplex2solver(b, descr.grid, temp_b);

	cg.Solve(temp_x, temp_b);

	simplexSupport::solver2simplex(descr.grid, temp_x, x);
}

void Poisson3D::solve_fast()
{
	//PoissonSmoothing *temp = new PoissonSmoothing((PoissonMapping*)cg.linear_mapping);
	for (int i = l - 2; i >= 0; i--) updateSubSystem(mg_descr[i], mg_descr[i + 1]);
	//if (cg.preconditioner) cg.preconditioner = nullptr;
	cg.preconditioner->update();

	cudaDeviceSynchronize();
	checkCudaErrors(cudaGetLastError());

	const PoissonDescriptor<3>& descr = mg_descr[l - 1];

	//simplexSupport::simplex2solver(b, descr.grid, temp_b);

	cg.Solve(temp_x, temp_b);

	//simplexSupport::solver2simplex(descr.grid, temp_x, x);
}

void Poisson3D::init_boundary(const FaceField<Scalar, 3>& face_vol, const Field<int, 3>& cell_fixed, bool _closed)
{
	PoissonDescriptor<3>& descr = mg_descr[l - 1];

	int size = descr.size;
	bool *fixed = new bool[size];
	Scalar *vol = new Scalar[descr.fsize];
	memset(fixed, 0, sizeof(bool)*descr.size);
	memset(vol, 0, sizeof(Scalar)*descr.fsize);

	simplexSupport::simplex2solver(cell_fixed, descr.grid, fixed, [=](int v)->bool { if (v) return true; else return false; });

	simplexSupport::simplex2solver(face_vol, descr.grid, vol);

	descr.setFixed(fixed);
	descr.setVol(vol);
	descr.toDevice();
	descr.finish();

	closed = _closed;

	delete[] fixed;
	delete[] vol;
}

void Poisson3D::update_b(const Field<Scalar, 3>& cell_b)
{
	b = cell_b;
}

