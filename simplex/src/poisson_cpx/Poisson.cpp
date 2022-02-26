#include "Poisson.h"

#include "PoissonMapping.h"
#include "PoissonMapping3D.h"
#include "PoissonLike.h"
#include "PoissonSubSystem.h"
#include "PoissonSubSystem3D.h"
#include "MultigridSolve.h"

#include "SimplexSupport.h"

#include "Timer.h"

template<int d>
void Poisson<d>::init(VectorDi _grid_size, Scalar _dx)
{
	grid_size = _grid_size;
	grid.Initialize(grid_size, _dx);
	int n = grid_size.minCoeff();
	if (d == 2) l = round(std::log2(n & -n)) - 2;
	else l = round(std::log2(n & -n)) - 1;

	x.Resize(grid_size);
	b.Resize(grid_size);

	PoissonDescriptor<d> descr = PoissonDescriptor<d>();
	descr.init(grid_size);

	temp_x = new Scalar[descr.size];
	temp_b = new Scalar[descr.size];

	cg = ConjugatedGradient();
	cg.CG_relative_residual_thres = 1e-3;
	//cg.linear_mapping = new PoissonDiagMapping(mapping);

	MultigridSolve *mg_solve = new MultigridSolve();
	mg_solve->init(l);

	mg_descr = new PoissonDescriptor<d>[l];
	mg_descr[l - 1] = descr;

	for (int i = l - 2; i >= 0; i--) mg_descr[i] = createSubSystem(mg_descr[i + 1]);

	using FixedMapping = typename If<d == 2, PoissonMappingFixed, PoissonMapping3DFixed >::Type;
	using DownSampler= typename If<d == 2, PoissonDownSample, PoissonDownSample3D >::Type;
	using UpSampler= typename If<d == 2, PoissonUpSample, PoissonUpSample3D >::Type;
	for (int i = 0; i < l; i++)
	{
		FixedMapping*t0 = new FixedMapping();
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

			DownSampler*t3 = new DownSampler(&mg_descr[i], &mg_descr[i - 1]);
			mg_solve->downSample[i - 1] = t3;

			UpSampler*t4 = new UpSampler(&mg_descr[i], &mg_descr[i - 1]);
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

	//cg.preconditioner = new DescentSmoother(new PoissonSmoothing((PoissonMapping*)cg.linear_mapping), 20, cg.linear_mapping);
	cg.preconditioner = mg_solve;
	cg.linear_mapping = mg_solve->mapping[l - 1];
	cg.Init(50);
}

template<int d>
void Poisson<d>::solve()
{
	for (int i = l - 2; i >= 0; i--) updateSubSystem(mg_descr[i], mg_descr[i + 1]);
	cg.preconditioner->update();

	cudaDeviceSynchronize();
	checkCudaErrors(cudaGetLastError());

	PoissonDescriptor<d>& descr = mg_descr[l - 1];

	simplexSupport::simplex2solver(b, descr.grid, temp_b);

	cg.Solve(temp_x, temp_b);

	simplexSupport::solver2simplex(descr.grid, temp_x, x);
}

template<int d>
void Poisson<d>::solve_fast()
{
	for (int i = l - 2; i >= 0; i--) updateSubSystem(mg_descr[i], mg_descr[i + 1]);
	cg.preconditioner->update();

	cudaDeviceSynchronize();
	checkCudaErrors(cudaGetLastError());

	PoissonDescriptor<d>& descr = mg_descr[l - 1];

	//simplexSupport::simplex2solver(b, descr.grid, temp_b);

	cg.Solve(temp_x, temp_b);

	//simplexSupport::solver2simplex(descr.grid, temp_x, x);
}

template<int d>
void Poisson<d>::init_boundary(const FaceField<Scalar, d>& face_vol, const Field<int, d>& cell_fixed, bool _closed)
{
	PoissonDescriptor<d>& descr = mg_descr[l - 1];

	int size = descr.size;
	bool *fixed = new bool[size];
	Scalar* vol = new Scalar[descr.fsize];
	memset(fixed, 0, sizeof(bool) * descr.size);
	memset(vol, 0, sizeof(Scalar) * descr.fsize);

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

template<int d>
void Poisson<d>::update_b(const Field<Scalar, d>& cell_b)
{
	b = cell_b;
}

template class Poisson<2>;
template class Poisson<3>;