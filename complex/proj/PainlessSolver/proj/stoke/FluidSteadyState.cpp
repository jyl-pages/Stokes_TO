#include "FluidSteadyState.h"
#include "StokeFlowMapping.h"
#include "StokeFlowSmoothing.h"
#include "StokeFlowDirectSolve.h"
#include "MultigridSolve.h"
#include "SimplexSupport.h"
#include "Timer.h"

void FluidSteadyState::init(int _grid_size)
{
	grid_size = _grid_size;
	l = round(std::log2(grid_size&-grid_size)) - 2;

	grid.Initialize(Vector2i(grid_size, grid_size), 1.);
	mac_grid.Initialize(grid);
	vel.Resize(Vector2i(grid_size, grid_size));
	p.Resize(Vector2i(grid_size, grid_size));

	MultigridSolve *mg_solve = new MultigridSolve();
	mg_solve->init(l);

	mg_descr = new StokeFlowDescriptor[l];
	StokeFlowDescriptor& descr = mg_descr[l - 1];

	descr.init(grid_size, grid_size);
	descr.lap_coeff = (Scalar)1;
	descr.toDevice();
	descr.finish();

	for (int i = l - 2; i >= 0; i--) mg_descr[i] = createSubSystem(mg_descr[i + 1]);

	for (int i = 0; i < l; i++)
	{
		StokeFlowMapping *t0 = new StokeFlowMapping();
		t0->init(&mg_descr[i]);
		mg_solve->mapping[i] = t0;

		if (i != 0)
		{

			StokeFlowSmoothing *t1 = new StokeFlowSmoothing();
			t1->init(t0);
			t1->order = 0;
			t1->iter_num = 1;
			mg_solve->preSmoother[i] = t1;
			//if (i != l - 2) mg_solve->preSmoother[i] = nullptr;

			StokeFlowSmoothing *t2 = new StokeFlowSmoothing();
			t2->init(t0);
			t2->order = 1;
			t2->iter_num = 1;
			mg_solve->postSmoother[i] = t2;
			//if (i != l - 2) mg_solve->postSmoother[i] = nullptr;

			StokeFlowDownSample *t3 = new StokeFlowDownSample();
			t3->init(&mg_descr[i], &mg_descr[i - 1]);
			mg_solve->downSample[i - 1] = t3;

			StokeFlowUpSample *t4 = new StokeFlowUpSample();
			t4->init(&mg_descr[i], &mg_descr[i - 1]);
			mg_solve->upSample[i - 1] = t4;
		}
	}

	StokeFlowDirectSolve *t5 = new StokeFlowDirectSolve();
	t5->init(&mg_descr[0]);
	t5->closed = closed;

	mg_solve->preSmoother[0] = t5;
	//mg_solve->preSmoother[0] = nullptr;
	mg_solve->postSmoother[0] = nullptr;

	mg_solve->finish();

	LinearMapping *outmost = mg_solve->mapping[l - 1];

	gmres.linear_mapping = outmost;
	gmres.preconditioner = mg_solve;
	//gmres.preconditioner = nullptr;
	//gmres.preconditioner = mg_solve->preSmoother[l - 1];
	gmres.Init(40);
	gmres.verbose = false;

	direct = t5;
}

void FluidSteadyState::init_boundary(const Field<Scalar, 2>& cell_vol, const Field<int, 2>& cell_fixed, \
	const FaceField<Scalar, 2>& face_vol, const FaceField<int, 2>& face_fixed, bool _closed)
{
	StokeFlowDescriptor& descr = mg_descr[l - 1];
	int size = descr.size;
	bool *fixed = new bool[size];
	Scalar *vol = new Scalar[size];
	memset(fixed, 0, sizeof(bool)*descr.size);
	memset(vol, 0, sizeof(Scalar)*descr.size);

	simplexSupport::simplex2solver(cell_vol, descr.grid, vol);
	simplexSupport::simplex2solver(cell_fixed, descr.grid, fixed, [=](int v)->bool { if (v) return true; else return false; });

	simplexSupport::simplex2solver(face_vol, descr.grid, vol + descr.grid.cell_size());
	simplexSupport::simplex2solver(face_fixed, descr.grid, fixed + descr.grid.cell_size(), [=](int v)->bool { if (v) return true; else return false; });

	descr.setFixed(fixed);
	descr.setVol(vol);

	closed = _closed;

	delete[] fixed;
	delete[] vol;
}

void FluidSteadyState::update_b(const Field<Scalar, 2>& cell_b, const FaceField<Scalar, 2>& face_b)
{
	StokeFlowDescriptor& descr = mg_descr[l - 1];
	int size = descr.size;
	Scalar *b = new Scalar[size];
	memset(b, 0, sizeof(Scalar)*descr.size);
	simplexSupport::simplex2solver(cell_b, descr.grid, b);
	simplexSupport::simplex2solver(face_b, descr.grid, b + descr.grid.cell_size());
	gmres.b.resize(size);
	memcpy(gmres.b.data(), b, sizeof(Scalar)*size);
	delete[] b;
}


void FluidSteadyState::update_penalty(const Field<Scalar, 2>& cell_penalty, const FaceField<Scalar, 2>& face_penalty)
{	
	StokeFlowDescriptor& descr = mg_descr[l - 1];

	int size = descr.size;
	Scalar *penalty= new Scalar[size];
	memset(penalty, 0, sizeof(Scalar)*descr.size);

	simplexSupport::simplex2solver(cell_penalty, descr.grid, penalty);
	simplexSupport::simplex2solver(face_penalty, descr.grid, penalty + descr.grid.cell_size());

	descr.setPenalty(penalty);
	descr.toDevice();
	descr.finish();
	for (int i = l - 2; i >= 0; i--) updateSubSystem(mg_descr[i], mg_descr[i + 1]);

	direct->finish();
	delete[] penalty;
}

void FluidSteadyState::solve()
{
	gmres.Solve();
	p.Fill((Scalar)0); vel.Fill((Scalar)0);

	const StokeFlowDescriptor& descr = mg_descr[l - 1];

	simplexSupport::solver2simplex(descr.grid, gmres.x.data(), p);
	simplexSupport::solver2simplex(descr.grid, gmres.x.data() + descr.grid.cell_size(), vel);

}
