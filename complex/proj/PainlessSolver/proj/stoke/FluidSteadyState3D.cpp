#include "FluidSteadyState3D.h"
#include "StokeFlowMapping3D.h"
//#include "StokeFlowSmoothing.h"
//#include "StokeFlowDirectSolve.h"
#include "StokeLike.h"
#include "StokeFlowSubSystem3D.h"
#include "MultigridSolve.h"
#include "SimplexSupport.h"
#include "LinearMappingDBG.h"
#include "Timer.h"

void FluidSteadyState3D::init(int _grid_size)
{
	grid_size = _grid_size;
	l = round(std::log2(grid_size&-grid_size)) - 1;

	Vector3i sz = Vector3i(grid_size, grid_size, grid_size);
	grid.Initialize(sz, 1.);
	mac_grid.Initialize(grid);
	vel.Resize(sz);
	p.Resize(sz);

	mg_descr = new StokeFlowDescriptor3D[l];
	StokeFlowDescriptor3D& descr = mg_descr[l - 1];

	descr.init(grid_size, grid_size, grid_size);
	descr.toDevice();
	descr.finish();

	MultigridSolve *mg_solve = new MultigridSolve();
	mg_solve->init(l);

	for (int i = l - 2; i >= 0; i--)
	{
		mg_descr[i] = createSubSystem(mg_descr[i + 1]);
		//mg_descr[i].toHost();
	}

	for (int i = 0; i < l; i++)
	{
		StokeFlowMapping3DFixed *t0 = new StokeFlowMapping3DFixed();
		t0->init(&mg_descr[i]);
		mg_solve->mapping[i] = t0;

		if (i != 0)
		{
			std::vector<LinearMapping*> pred = StokeLike::GenerateJacobiPreconditioner(mg_descr[i].grid, t0, 2);
			mg_solve->preSmoother[i] = pred[0];
			//if (i != 0) mg_solve->preSmoother[i] = nullptr;

			//mg_solve->postSmoother[i] = mg_solve->preSmoother[i];
			mg_solve->postSmoother[i] = pred[1];
			//if (i != 0) mg_solve->postSmoother[i] = nullptr;

			StokeFlowDownSample3D *t3 = new StokeFlowDownSample3D();
			t3->init(&mg_descr[i], &mg_descr[i - 1]);
			mg_solve->downSample[i - 1] = t3;

			StokeFlowUpSample3D *t4 = new StokeFlowUpSample3D();
			t4->init(&mg_descr[i], &mg_descr[i - 1]);
			mg_solve->upSample[i - 1] = t4;
		}

		else
		{
			mg_solve->preSmoother[i] = StokeLike::GenerateDirectSolve(mg_descr[i].grid, t0, closed);
			//mg_solve->preSmoother[i] = nullptr;
			mg_solve->postSmoother[i] = nullptr;
		}
	}

	//StokeFlowDirectSolve *t5 = new StokeFlowDirectSolve();
	//t5->init(&mg_descr[0]);
	//t5->closed = closed;

	//mg_solve->preSmoother[0] = t5;
	////mg_solve->preSmoother[0] = nullptr;
	//mg_solve->postSmoother[0] = nullptr;

	mg_solve->finish();

	gmres.linear_mapping = mg_solve->mapping[l - 1];
	//gmres.linear_mapping = StokeLike::GenerateJacobiPreconditioner(descr.grid, outmost);
	gmres.preconditioner = mg_solve;
	//gmres.preconditioner = nullptr;
	//gmres.preconditioner = mg_solve->preSmoother[l - 1];
	gmres.Init(50);

	//direct = t5;
}

void FluidSteadyState3D::init_boundary(const Field<Scalar, 3>& cell_vol, const Field<int, 3>& cell_fixed, \
	const FaceField<Scalar, 3>& face_vol, const FaceField<int, 3>& face_fixed, bool _closed)
{
	StokeFlowDescriptor3D& descr = mg_descr[l - 1];
	bool *fixed = new bool[descr.size];
	Scalar *vol = new Scalar[descr.tsize];
	memset(fixed, 0, sizeof(bool)*descr.size);
	memset(vol, 0, sizeof(Scalar)*descr.tsize);

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

void FluidSteadyState3D::update_b(const Field<Scalar, 3>& cell_b, const FaceField<Scalar, 3>& face_b)
{
	StokeFlowDescriptor3D& descr = mg_descr[l - 1];
	int size = descr.size;
	Scalar *b = new Scalar[size];
	memset(b, 0, sizeof(Scalar)*descr.size);
	simplexSupport::simplex2solver(cell_b, descr.grid, b);
	simplexSupport::simplex2solver(face_b, descr.grid, b + descr.grid.cell_size());
	gmres.b.resize(size);
	memcpy(gmres.b.data(), b, sizeof(Scalar)*size);
	delete[] b;
}


void FluidSteadyState3D::update_penalty(const Field<Scalar, 3>& cell_penalty, const FaceField<Scalar, 3>& face_penalty)
{
	StokeFlowDescriptor3D& descr = mg_descr[l - 1];

	int size = descr.size;
	Scalar *penalty = new Scalar[size];
	memset(penalty, 0, sizeof(Scalar)*descr.size);

	simplexSupport::simplex2solver(cell_penalty, descr.grid, penalty);
	simplexSupport::simplex2solver(face_penalty, descr.grid, penalty + descr.grid.cell_size());

	descr.setPenalty(penalty);
	descr.toDevice();
	descr.finish();

	//for (int i = l - 2; i >= 0; i--) updateSubSystem(mg_descr[i], mg_descr[i + 1]);

	//direct->finish();

	delete[] penalty;
}

void FluidSteadyState3D::solve()
{
	StokeFlowDescriptor3D& descr = mg_descr[l - 1];
	for (int i = l - 2; i >= 0; i--) updateSubSystem(mg_descr[i], mg_descr[i + 1]);
	gmres.preconditioner->update();
	//LinearMappingDBG::LinearMapping2SparseMatrix("matrix.txt", gmres.linear_mapping, 1e-3);

	gmres.Solve();
	p.Fill((Scalar)0); vel.Fill((Scalar)0);


	//LinearMappingDBG::LinearMapping2DenseMatrix("matrix.txt", gmres.linear_mapping);
	//LinearMappingDBG::LinearMapping2SparseMatrix("matrix.txt", gmres.preconditioner, 1e-2);
	//LinearMappingDBG::LinearMapping2SparseMatrix("matrix.txt", mg_solve->mapping[0], 1e-2);

	simplexSupport::solver2simplex(descr.grid, gmres.x.data(), p);
	simplexSupport::solver2simplex(descr.grid, gmres.x.data() + descr.grid.cell_size(), vel);

}
