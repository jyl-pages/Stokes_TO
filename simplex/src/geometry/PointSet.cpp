//////////////////////////////////////////////////////////////////////////
// Point set
// Copyright (c) (2018-), Bo Zhu, Hui Wang, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "PointSet.h"

template<int d>
void PointSet<d>::Update_Nbs()
{
	//update data structure nbs_searcher
	std::function<bool(const int)> valid = [&](const int idx)->bool {return points->I(idx) != -1; };
	nbs_searcher->Update_Points(points->XRef(), valid);
	//use updated nbs_searcher to update tang_nbs
	int pn = points->Size();
	tang_nbs.resize(pn);
#pragma omp parallel for
	for (int i = 0; i < pn; i++) {
		if (!Valid_Particle(i)) continue;
		const VectorD &pos = points->X(i);
		const MatrixD& local_frame = points->E(i);
		std::function<bool(const int)> filter_func = [&](const int p) {return Is_Tangential_Neighbor(pos, local_frame, p); };
		nbs_searcher->Find_Neighbors(pos, v_r, filter_func, tang_nbs[i]);
	}
}

template<int d>
bool PointSet<d>::Is_Tangential_Neighbor(const VectorD& pos, const MatrixD& local_frame, const int p) const
{
	////check angle
	VectorD dir = local_frame.col(d - 1);
	VectorD dir_p = points->E(p).col(d - 1);
	real dot = dir.dot(dir_p);
	if (dot < t_dot) return false;	////skip the points with large angles

	////check distance
	VectorD u = points->X(p) - pos;
	VectorT t = Project_To_TPlane(u, local_frame);
	return t.norm() < t_r;
}

template<int d>
ArraySlice<int> PointSet<d>::Record_Tangential_Nbs(const VectorD& pos, const MatrixD& local_frame) const
{
	std::function<bool(const int)> filter_func = [&](const int p) {return Is_Tangential_Neighbor(pos, local_frame, p); };
	return nbs_searcher->Record_Neighbors(pos, v_r, filter_func);
}

template<int d>
Array<int> PointSet<d>::Find_Tangential_Nbs(const VectorD& pos, const MatrixD& local_frame) const
{
	std::function<bool(const int)> filter_func = [&](const int p) {return Is_Tangential_Neighbor(pos, local_frame, p); };
	return nbs_searcher->Find_Neighbors(pos, v_r, filter_func);
}

template<int d>
void PointSet<d>::Update_Local_Frame(const real dt)
{
	int p_n = points->Size();
	for (int i = 0; i < p_n; i++) {
		if (!Valid_Particle(i))continue;

		if constexpr (d == 2) {
			VectorD e0 = points->E(i).col(0);
			VectorD e1 = points->E(i).col(1);

			////return u \dot e1 in tangential space
			std::function<real(const int)> V_Normal_Val = [&](const int idx)->real
			{return points->V(idx).dot(e1); };

			VectorT dvdx = Grad_TPlane(i, V_Normal_Val);
			std::cout << dvdx[0] << ", ";
			e0 = Rotation2(dvdx[0] * dt) * e0;
			e1 = Rotation2(dvdx[0] * dt) * e1;
		}
		else if constexpr (d == 3) {
			////TOIMPROVE: use rotation matrix or quaternion
			VectorD e0 = points->E(i).col(0);
			VectorD e1 = points->E(i).col(1);
			VectorD e2 = points->E(i).col(2);

			std::function<real(const int)> U0_Val = [&](const int idx)->real
			{return points->V(idx).dot(points->E(i).col(0)); };
			VectorT du0 = Grad_TPlane(i, U0_Val);

			std::function<real(const int)> U1_Val = [&](const int idx)->real
			{return points->V(idx).dot(points->E(i).col(1)); };
			VectorT du1 = Grad_TPlane(i, U1_Val);

			std::function<real(const int)> U2_Val = [&](const int idx)->real
			{return points->V(idx).dot(points->E(i).col(2)); };
			VectorT du2 = Grad_TPlane(i, U2_Val);

			VectorD de0 = (real).5 * (du1[0] - du0[1]) * e1 + (du2[0]) * e2;
			VectorD de1 = (real).5 * (du0[1] - du1[0]) * e0 + (du2[1]) * e2;
			VectorD de2 = (-du2[0]) * e0 + (-du2[1]) * e1;

			e0 += de0 * dt; e0.normalize();
			e1 += de1 * dt; e1.normalize();
			e2 += de2 * dt; e2.normalize();

			points->E(i).col(0) = e0;
			points->E(i).col(1) = e1;
			points->E(i).col(2) = e2;
		}
	}
}

template<int d>
void PointSet<d>::Initialize_Local_Frame(const VectorD& v, const Array<int>& nbs,MatrixD& local_frame){

	////calculate the PCA normal
	VectorD xp = VectorD::Zero();
	real w = (real)0;
	for (int i = 0; i < (int)nbs.size(); i++) {
		////here we use volumetric distance instead of tangential distance
		real dis = (points->X(i) - v).norm();
		real w0 = W_PCA(dis);
		xp += w0 * points->X(i);
		w += w0;
	}
	if (w != (real)0)xp /= w;
	MatrixD C = MatrixD::Zero();
	real wc = (real)0;
	for (int i = 0; i < (int)nbs.size(); i++) {
		const VectorD& xi = points->X(i);
		real dis = (xi - xp).norm();
		real w0 = W_PCA(dis);
		C += w0 * (xi - xp) * (xi - xp).transpose();
		wc += w0;
	}
	if (wc != (real)0)C /= wc;
	VectorD normal = AuxFunc::Min_Eigenvector(C);
	// if (normal.dot(points->Normal(p)) < (real)0)normal *= (real)-1;
	if (normal.dot(Normal(v)) < (real)0)normal *= (real)-1; ////TOTEST

	////update local frame according to the PCA normal
	if constexpr (d == 2) {
		VectorD tang = -AuxFunc::Orthogonal_Vector(normal);
		local_frame.col(0) = tang.normalized();
		local_frame.col(1) = normal.normalized();
	}
	else if constexpr (d == 3) {
		VectorD t1 = -AuxFunc::Orthogonal_Vector(normal);
		VectorD t2 = t1.cross(normal);
		local_frame.col(0) = t1.normalized();
		local_frame.col(1) = t2.normalized();
		local_frame.col(2) = normal.normalized();
	}
}

template<int d>
void PointSet<d>::Reinitialize_Local_Frames()
{
	////update local normals with PCA
	int pn = points->Size();
#pragma omp parallel for
	for (int p = 0; p < pn; p++) {
		if (!Valid_Particle(p))continue;

		////calculate the PCA normal
		size_t nbs_num = tang_nbs[p].size();
		VectorD xp = VectorD::Zero();
		real w = (real)0;
		for (size_t i = 0; i < nbs_num; i++) {
			int q = tang_nbs[p][i];
			////here we use volumetric distance instead of tangential distance
			real dis = (points->X(q) - points->X(p)).norm();
			real w0 = W_PCA(dis);
			xp += w0 * points->X(q);
			w += w0;
		}
		if (w != (real)0)xp /= w;
		MatrixD C = MatrixD::Zero();
		real wc = (real)0;
		for (int i = 0; i < nbs_num; i++) {
			int q = tang_nbs[p][i];
			const VectorD& xq = points->X(q);
			real dis = (xq - xp).norm();
			real w0 = W_PCA(dis);
			C += w0 * (xq - xp) * (xq - xp).transpose();
			wc += w0;
		}
		if (wc != (real)0)C /= wc;
		VectorD normal = AuxFunc::Min_Eigenvector(C);
		if (normal.dot(points->Normal(p)) < (real)0)normal *= (real)-1;

		////update local frame according to the PCA normal
		if constexpr (d == 2) {
			VectorD tang = -AuxFunc::Orthogonal_Vector(normal);
			points->E(p).col(0) = tang.normalized();
			points->E(p).col(1) = normal.normalized();
		}
		else if constexpr (d == 3) {
			VectorD t1 = -AuxFunc::Orthogonal_Vector(normal);
			VectorD t2 = t1.cross(normal);
			points->E(p).col(0) = t1.normalized();
			points->E(p).col(1) = t2.normalized();
			points->E(p).col(2) = normal.normalized();
		}
	}

		//#pragma omp parallel for
		//for (int p = 0; p < pn; p++) {
		//	
		//}

	////Correct local frame to align with the tangent space
	const bool use_local_frame_correction = true;
	if(use_local_frame_correction) {////
		#pragma omp parallel for
		for (int p = 0; p < pn; p++) {
			if (!Valid_Particle(p))continue;
			////return h in tangential space
			std::function<real(const int)> H_Val = [=](const int idx)->real {
				VectorD u = points->X(idx) - points->X(p); real h = Project_To_TPlane_H(u, points->E(p)); return h; };

			VectorT dzdx = Grad_TPlane(p, H_Val);

			if constexpr (d == 2) {
				VectorD tang = (points->E(p) * VectorD(1, dzdx(0))).normalized();
				VectorD normal = -AuxFunc::Orthogonal_Vector(tang);

				points->E(p).col(0) = tang.normalized();
				points->E(p).col(1) = normal.normalized();
			}
			else if constexpr (d == 3) {
				VectorD t0 = points->E(p) * VectorD(1, 0, dzdx(0));
				VectorD t1 = points->E(p) * VectorD(0, 1, dzdx(1));
				VectorD normal = t1.cross(t0).normalized();

				t0 = -AuxFunc::Orthogonal_Vector(normal);
				t1 = t0.cross(normal);
				points->E(p).col(0) = t0.normalized();
				points->E(p).col(1) = t1.normalized();
				points->E(p).col(2) = normal.normalized();
			}
		}
	}
}

template<int d>
Vector<real,d> PointSet<d>::Normal(const VectorD& pos) const
{
	int closest_p = Closest_Point(pos);
	Array<int> nbs = Find_Tangential_Nbs(pos, points->E(closest_p));
	size_t nb_n = nbs.size();
	if (nb_n == 0) {
		return points->Normal(closest_p);
	}

	VectorD nml = VectorD::Zero();
	for (int i = 0; i < nb_n; i++) {
		int p = nbs[i];
		real dis = (pos - points->X(p)).norm();
		real w0 = W_PCA(dis);
		nml += w0 * Normal(p);
	}
	return nml.normalized();
}

template<int d>
Matrix<real, d> PointSet<d>::Local_Frame(const VectorD& pos) const
{
	VectorD normal = Normal(pos);
	MatrixD e;
	if constexpr (d == 2) {
		VectorD tang = -AuxFunc::Orthogonal_Vector(normal);
		e.col(0) = tang.normalized();
		e.col(1) = normal.normalized();
	}
	else if constexpr (d == 3) {
		VectorD t1 = -AuxFunc::Orthogonal_Vector(normal);
		VectorD t2 = t1.cross(normal);
		e.col(0) = t1.normalized();
		e.col(1) = t2.normalized();
		e.col(2) = normal.normalized();
	}
	return e;
}

template<int d>
void PointSet<d>::Update_Metric_Tensor()
{
	int p_n = points->Size();
	for (int i = 0; i < p_n; i++) {
		if (!Valid_Particle(i))continue;

		////return h in tangential space
		std::function<real(const int)> H_Val = [=](const int idx)->real
		{VectorD u = points->X(idx) - points->X(i); real h = Project_To_TPlane_H(u, points->E(i)); return h; };

		VectorT dzdx = Grad_TPlane(i, H_Val);
		points->dH(i) = dzdx;
		points->G(i) = Metric_Tensor(dzdx);
	}
}

template<int d>
Matrix<real, 1> PointSet<d>::Metric_Tensor(const Vector1& dzdx) const
{
	Matrix<real, 1> mt; mt << 1 + pow(dzdx[0], 2);
	return mt;
}

template<int d>
Matrix<real, 2> PointSet<d>::Metric_Tensor(const Vector2& dzdx) const
{
	Matrix<real, 2> mt;
	real g11 = 1 + pow(dzdx[0], 2);
	real g22 = 1 + pow(dzdx[1], 2);
	real g12 = dzdx[0] * dzdx[1];
	mt << g11, g12, g12, g22;
	return mt;
}

template<int d>
Vector<real, d> PointSet<d>::Local_Coords(const VectorD& u, const MatrixD& e)
{
	//[u.dot(e.col(0)), u.dot(e.col(1)), ...]
	return u.transpose() * e;
}

template<int d>
void PointSet<d>::Local_Coords(const VectorD& u, const MatrixD& e, VectorT& t, real& z)
{
	VectorD local_coords = Local_Coords(u, e);
	t = AuxFunc::V<d - 1>(local_coords);
	z = local_coords[d - 1];
}

template<int d>
Vector<real, d - 1> PointSet<d>::Project_To_TPlane(const VectorD& u, const MatrixD& e)
{
	VectorT t_coords;
	for (int i = 0; i < d - 1; i++) {
		t_coords[i] = u.dot(e.col(i));
	}
	return t_coords;
}

template<int d>
real PointSet<d>::Project_To_TPlane_H(const VectorD& u, const MatrixD& e)
{
	return u.dot(e.col(d - 1));
}

template<int d>
int PointSet<d>::Nearest_Geometry(VectorD& pos, MatrixD& local_frame, int minimum_eqns, bool verbose) const
{
	////The current implementation only conducts one iteration. The function can be called for multiple times for an iterative projection.
	using namespace LeastSquares;
	local_frame = Local_Frame(pos);
	Array<int> nbs = Find_Tangential_Nbs(pos, local_frame);
	MLS<d - 1, 2> ls;
	size_t n = nbs.size(); Array<real> data(n * (size_t)d);
	if (n < minimum_eqns) {
		int closest_p = Closest_Point(pos);
		VectorD norm = points->Normal(closest_p);
		real offset = norm.normalized().dot(pos - points->X(closest_p));
		pos = pos - offset * norm;
		if(verbose)	std::cout << "[Warning]PointSet<d>::Project_To_Surface: " << pos.transpose() << " find only " << n << " neighbors against minimum " << minimum_eqns << "\n";
		return 1;
	}
	else {
		for (size_t i = 0; i < n; i++) {
			int p = nbs[i];
			VectorD u = points->X(p) - pos;
			VectorD th = Local_Coords(u, local_frame);
			for (size_t j = 0; j < d; j++)data[i * (size_t)d + j] = th[j];
		}
		ls.Fit(data.data(), (int)n, 0, 0);
		//pos = pos + local_frame.col(d - 1) * ls(VectorT::Zero());
		pos = pos + local_frame.col(d - 1) * ls(VectorT::Zero());
		return 0;
	}
}

template<int d>
Vector<real, d> PointSet<d>::Project_To_Surface(const VectorD& pos) const
{
	VectorD proj_pos = pos;
	MatrixD frame;
	Nearest_Geometry(proj_pos, frame);
	return proj_pos;
}

template<int d>
real PointSet<d>::Unsigned_Distance(const VectorD& pos) const
{
	VectorD proj_pos = Project_To_Surface(pos);
	VectorD v = proj_pos - pos;
	return v.norm();
}

template<int d>
real PointSet<d>::Signed_Distance(const VectorD& pos) const
{
	VectorD proj_pos = Project_To_Surface(pos);
	VectorD v = proj_pos - pos;
	VectorD normal = Normal(proj_pos);
	real sign = v.normalized().dot(normal) > (real)0 ? (real)1 : (real)-1;
	return sign * v.norm();
}

template<int d>
real PointSet<d>::Truncated_Unsigned_Distance(const VectorD& pos, const real& max_dist, int minimum_eqns) const
{
	int closest_p = Closest_Point(pos);
	real closest_dist = (points->X(closest_p) - pos).norm();
	if (max_dist > 0 && closest_dist > max_dist && closest_dist > v_r) {
		return max_dist;
	}
	VectorD proj_pos = pos; MatrixD frame;
	Nearest_Geometry(proj_pos, frame, minimum_eqns);
	real proj_dist = (proj_pos - pos).norm();
	if (max_dist > 0 && proj_dist > max_dist) {
		return max_dist;
	}
	return proj_dist;
}

template<int d>
real PointSet<d>::Truncated_Signed_Distance(const VectorD& pos, const real& max_dist, int minimum_eqns) const
{
	real phi = Truncated_Unsigned_Distance(pos, max_dist, minimum_eqns);
	VectorD proj_pos = pos; MatrixD frame;
	Nearest_Geometry(proj_pos, frame, minimum_eqns);
	VectorD v = proj_pos - pos;
	VectorD normal = frame.col(d - 1);
	real sign = v.normalized().dot(normal) > (real)0 ? (real)1 : (real)-1;
	return sign * phi;
}

template<int d>
Vector<real,d> PointSet<d>::Curvature_Vector(const int i)
{
	std::function<real(const int)> Point_Position[d];

	Point_Position[0] = [&](const int idx)->real
	{return points->X(idx)[0]; };
	Point_Position[1] = [&](const int idx)->real
	{return points->X(idx)[1]; };
	if constexpr (d == 3)
		Point_Position[2] = [&](const int idx)->real
	{return points->X(idx)[2]; };

	VectorD curv_vec;
	for (int a = 0; a < d; a++) {
		curv_vec[a] = Laplacian(i, Point_Position[a]);
	}
	return curv_vec;

	// VectorX geo_coeff=Fit_Local_Geometry_MLS(i);
	// MatrixT g;
	// if constexpr(d==2){
	// 	g(0,0)=geo_coeff(1)*geo_coeff(1)+1;
	// }else if constexpr(d==3){
	// 	g(0,0)=geo_coeff(1)*geo_coeff(1)+1;
	// 	g(0,1)=geo_coeff(1)*geo_coeff(2);
	// 	g(1,0)=geo_coeff(1)*geo_coeff(2);
	// 	g(1,1)=geo_coeff(2)*geo_coeff(2)+1;}
	// MatrixT g_inv=g.inverse();
	// // std::cout<<g_inv<<std::endl;

	// std::function<real(const int)> Point_Position[d];
	// VectorD curv=VectorD::Zero();
	// for(int a=0;a<d;a++){
	// 	Point_Position[a]=[&](const int idx)->real{return this->points->X(idx)[a];};
	// 	VectorX c=this->Fit_Local_MLS(i,Point_Position[a]);
	// 	// curv[a]=2*g_inv(0,0)*c(3)+2*g_inv(1,1)*c(4)+(g_inv(0,1)+g_inv(1,0))*c(5);
	// 	curv[a]=2*c(3)+2*c(4);
	// }
	// return curv;
}

template<int d>
Vector<real, d> PointSet<d>::Curvature_Vector(const VectorD& pos, const MatrixD& lf)
{
	std::function<real(const int)> Point_Position[d];

	Point_Position[0] = [&](const int idx)->real
	{return points->X(idx)[0]; };
	Point_Position[1] = [&](const int idx)->real
	{return points->X(idx)[1]; };
	if constexpr (d == 3)
		Point_Position[2] = [&](const int idx)->real
	{return points->X(idx)[2]; };

	VectorD curv_vec;
	for (int a = 0; a < d; a++) {
		curv_vec[a] = Laplacian(pos, lf, Point_Position[a]);
	}
	return curv_vec;
}

template<int d>
Vector<real, d - 1> PointSet<d>::Local_Dzdx(const int i, const int j)
{
	if constexpr (d == 2) {
		const VectorD ti = points->E(i).col(0);
		const VectorD ni = points->E(i).col(1);
		const VectorD nj = points->E(j).col(2);
		const real nj_ti = nj.dot(ti);
		real nj_ni = nj.dot(ni);

		if (abs(nj_ni) < 1e-8) {
			if (nj_ni < 0) nj_ni = 1e-8;
			else nj_ni = -1e-8;
		}
		return VectorT(nj_ti / nj_ni);
	}
	else if constexpr (d == 3) {
		const VectorD ti0 = points->E(i).col(0);
		const VectorD ti1 = points->E(i).col(1);
		const VectorD ni = points->E(i).col(2);
		const VectorD nj = points->E(j).col(2);
		const real nj_ti0 = nj.dot(ti0);
		const real nj_ti1 = nj.dot(ti1);
		real nj_ni = nj.dot(ni);

		if (abs(nj_ni) < 1e-8) {
			if (nj_ni < 0) nj_ni = 1e-8;
			else nj_ni = -1e-8;
		}
		return VectorT(nj_ti0 / nj_ni, nj_ti1 / nj_ni);
	}
}

template<int d>
Vector<real, d - 1> PointSet<d>::Local_Dzdx(const MatrixD& lf, const int j)
{
	if constexpr (d == 2) {
		const VectorD ti = lf.col(0);
		const VectorD ni = lf.col(1);
		const VectorD nj = points->E(j).col(2);
		const real nj_ti = nj.dot(ti);
		real nj_ni = nj.dot(ni);

		if (abs(nj_ni) < 1e-8) {
			if (nj_ni < 0) nj_ni = 1e-8;
			else nj_ni = -1e-8;}
		return VectorT(nj_ti / nj_ni);
	}
	else if constexpr (d == 3) {
		const VectorD ti0 = lf.col(0);
		const VectorD ti1 = lf.col(1);
		const VectorD ni = lf.col(2);
		const VectorD nj = points->E(j).col(2);
		const real nj_ti0 = nj.dot(ti0);
		const real nj_ti1 = nj.dot(ti1);
		real nj_ni = nj.dot(ni);

		if (abs(nj_ni) < 1e-8) {
			if (nj_ni < 0) nj_ni = 1e-8;
			else nj_ni = -1e-8;}
		return VectorT(nj_ti0 / nj_ni, nj_ti1 / nj_ni);
	}
}

template<int d>
real PointSet<d>::Calculate_Number_Density(const Array<VectorD>& X, Array<real>& nden)
{
	int pn = points->Size(); real avg = (real)0; int n = 0; nden.resize(pn); //n is the number of active particles
#pragma omp parallel for
	for (int i = 0; i < pn; i++) {
		if (!Valid_Particle(i))continue;
		real nd = (real)0;
		size_t nb_n = tang_nbs[i].size();
		for (size_t k = 0; k < nb_n; k++) {
			int j = tang_nbs[i][k];
			VectorT lr_ij = Project_To_TPlane(X[i] - X[j], points->E(i));
			nd += t_kernel.Weight(d, lr_ij.norm(), KernelType::POLY6);
		}
		nden[i] = nd; avg += nd; n++;
	}
	if (n > 0)avg /= (real)n; return avg;
}

template<int d>
real PointSet<d>::Calculate_Number_Density(const Array<VectorD>& X, const VectorD& pos) const
{
	MatrixD local_frame = Local_Frame(pos);
	Array<int> nbs_arr = Find_Tangential_Nbs(pos, local_frame);
	real nd = (real)0; size_t nb_n = nbs_arr.size();
	for (size_t k = 0; k < nb_n; k++) {
		size_t j = nbs_arr[k];
		VectorT lr_ij = Project_To_TPlane(pos - X[j], local_frame);
		nd += t_kernel.Weight(d, lr_ij.norm(), KernelType::POLY6);
	}
	return nd;
}

template<int d>
void PointSet<d>::Point_Reseeding()
{
	real reseeding_nden = (real)0.6 * init_avg_nden;
	real deleting_nden = (real)2 * init_avg_nden;
	int pn = points->Size();

	////find all reseeding pairs
	Array<Vector2i> reseeding_pairs;
	//Array<int> deleting_idx;
	for (int i = 0; i < pn; i++) {
		if (!Valid_Particle(i))continue;
		size_t nb_n = tang_nbs[i].size();
		for (size_t k = 0; k < nb_n; k++) {
			int j = tang_nbs[i][k];
			if (j <= i)continue;//avoid counting itself
			VectorD mid_pos = (real).5 * (points->X(i) + points->X(j));
			real mid_nden = Calculate_Number_Density(points->XRef(), mid_pos);
			if (mid_nden < reseeding_nden) {
				reseeding_pairs.push_back(Vector2i(i, j));
			}
		}
		real i_nden = Calculate_Number_Density(points->XRef(), points->X(i));
		//if(i_nden>deleting_nden)
		//	deleting_idx.push_back(i);
	}

	////reseed points according to the selected pairs
	reseeded_points.clear();
	int new_size = pn + (int)reseeding_pairs.size();
	points->Resize(new_size);
	for (int k = 0; k < reseeding_pairs.size(); k++) {
		int i = reseeding_pairs[k][0]; int j = reseeding_pairs[k][1];
		//int p=Add_Particle(); ////some bug with Add_Particles()
		int p = pn + k;
		reseeded_points.push_back(p);
		points->X(p) = (real).5 * (points->X(i) + points->X(j));
		points->V(p) = (real).5 * (points->V(i) + points->V(j));
		points->E(p) = Local_Frame(points->X(p));
		points->M(p) = (real)1;
		points->I(p) = 0;
	}
	if (reseeding_pairs.size() > 0) Update();
	////deleting points with high density
	//for (int k=0; k<deleting_idx.size(); k++)
	//{Remove_Particle(k);
	//std::cout<<"particle removed"<<std::endl;
	//}
}

template<int d>
void PointSet<d>::Point_Relaxation()
{
	int pn = points->Size();
	Array<real> nden(pn);
	Array<VectorT> local_f(pn);
	Array<VectorD> relaxed_pos(pn);

	real kp = (real)1e3;
	real nden_0 = (real).8 * init_avg_nden;	////set the nden_0 to be smaller than the initial average to force points to push each other
	real one_over_m = (real)1;
	real vis = (real)1;
	real dt = (real).02;
	int relax_iter_num = 4;

	auto P = [&](const int idx)->real {return kp * (nden[idx] / nden_0 - (real)1.); };	////lambda function to access P
	auto Vol = [&](const int idx)->real {return (real)1 / nden[idx]; };
	relaxed_pos = points->XRef();	////initialize relaxed_pos by *copying* data from points

	for (int iter = 0; iter < relax_iter_num; iter++) {
		////update number density
		for (int i = 0; i < pn; i++) {
			if (points->I(i) == -1)continue;
			real nd = (real)0;
			int nb_n = (int)tang_nbs[i].size();
			for (int k = 0; k < nb_n; k++) {
				int j = tang_nbs[i][k];
				VectorT lr_ij = Project_To_TPlane(relaxed_pos[i] - relaxed_pos[j], points->E(i));
				nd += t_kernel.Weight(d, lr_ij.norm(), KernelType::POLY6);
			}
			nden[i] = nd;
		}

		////update forces
		for (int i = 0; i < pn; i++) {
			if (points->I(i) == -1)continue;
			VectorT lf = VectorT::Zero();
			size_t nb_n = tang_nbs[i].size();
			for (size_t k = 0; k < nb_n; k++) {
				int j = tang_nbs[i][k];
				VectorT lr_ij = Project_To_TPlane(relaxed_pos[i] - relaxed_pos[j], points->E(i));
				real lr2 = lr_ij.squaredNorm();
				VectorT lf_p = -(P(i) * pow(Vol(i), 2) + P(j) * pow(Vol(j), 2)) * t_kernel.Grad<d - 1>(lr_ij, KernelType::SPIKY);
				lf += lf_p;
			}
			local_f[i] = one_over_m * lf;
		}

		////time integration
		for (int i = 0; i < pn; i++) {
			if (points->I(i) == -1)continue;
			VectorD delta_x; Unproject_To_World(VectorT(local_f[i] * dt * dt), points->E(i), delta_x);
			relaxed_pos[i] += delta_x;
			relaxed_pos[i] = Project_To_Surface(relaxed_pos[i]);
		}
	}

	////update positions
	for (int i = 0; i < pn; i++) {
		if (points->I(i) == -1)continue;
		points->X(i) = relaxed_pos[i];
	}
}

template<int d>
real PointSet<d>::Ray_Casting_MLS(VectorD ray_p, VectorD ray_d, real max_distance, MatrixD& lf_o, VectorX& c_o)
{
	VectorD tmp_p = ray_p;
	real t = (real)0;
	int cp_idx = -1;
	lf_o = MatrixD::Zero();
	c_o = VectorX(6);

	while (t < max_distance) {

		//// 1. find the closest point
		while (t < max_distance) {
			cp_idx = Closest_Point(tmp_p);
			if (cp_idx < 0) {//// if find no cp, move v_r forward
				t += v_r;
			} //// WH: I guess v_r<radius in cp searching?
			else {//// if find cp, move d(cp)-v_r forward
				VectorD cp = points->X(cp_idx);
				real step = (tmp_p - cp).norm() - v_r;
				if (step <= 1e-8) break;
				t += step;
			}
			tmp_p = ray_p + ray_d * t;
		}
		//std::cout<<"- Find cp "<<cp_idx<<" ("<<points->X(cp_idx).transpose()<<"), t="<<t<<", tmp_p ("<<tmp_p.transpose()<<")"<<std::endl;
		if (t > max_distance) return max_distance;

		//// 2. compute ray intersection and project
		real et = (real)0;
		VectorD hit_p = tmp_p; //// hit_p=ray_p+ray_d*(t+et)
		VectorD proj_p = points->X(cp_idx); //// proj_p=Project(hit_p)
		//// iteratively compute the ray intersection and projection
		while (true) {
			//// compute intersection
			{
				MatrixD lf; VectorX c = Fit_Local_Geometry_MLS(proj_p, lf);
				VectorD lp = Local_Coords(tmp_p - proj_p, lf);
				VectorD ld = Local_Coords(ray_d, lf);
				real new_et = Solve_Intersection(lp, ld, c);
				if (new_et < -v_r) break;
				et = new_et;
				hit_p = ray_p + ray_d * (t + et);
			}
			// std::cout<<"- hit_p=("<<hit_p.transpose()<<")"<<std::endl;
			//// project
			{
				MatrixD lf; VectorX c = Fit_Local_Geometry_MLS(hit_p, lf);
				if (abs(c(0)) < v_r * 0.1f) { c_o = c; lf_o = lf; return (real)t + et; }
				if (abs(c(0)) > v_r) { et = 0; break; }
				proj_p = hit_p - c(0) * lf.col(2);
			}
			// std::cout<<"- proj_p=("<<proj_p.transpose()<<")"<<std::endl;
		}

		//// update t and tmp_p
		if (et < v_r * 0.1) et = v_r * 0.1;
		t = t + et; tmp_p = ray_p + ray_d * t;
	}
	return max_distance;
}

template<int d>
void PointSet<d>::Solve_Ray_Sphere_Intersection(const VectorD& p, const VectorD& dir, const VectorD& s, const real r, real& norm_dist, real& tang_dist, real& delta_dist)
{
	VectorD rel = s - p;
	tang_dist = rel.dot(dir);
	VectorD norm = rel - tang_dist * dir;//// output
	norm_dist = norm.norm();
	real dist_2 = r * r - norm.squaredNorm();
	delta_dist = dist_2 < 0 ? -1 : sqrt(dist_2);
}

template<int d>
real PointSet<d>::Solve_Intersection(const VectorD& p, const VectorD& dir, const VectorX& coeff)
{
	if constexpr (d == 2)
		std::cerr << "Error: [PointSet] Unimplemented template<d=2>Solve_Intersection(...)" << std::endl;
	if constexpr (d == 3) {
		//// construct at^2+bt+c=0
		real a = coeff[3] * dir[0] * dir[0] + coeff[4] * dir[1] * dir[1] + coeff[5] * dir[0] * dir[1];//// a=(b3d0d0+b5d1d1+b4d0d1)
		real b = coeff[1] * dir[0] + coeff[2] * dir[1] + 2 * coeff[3] * p[0] * dir[0] + 2 * coeff[4] * p[1] * dir[1] + coeff[5] * p[0] * dir[1] + coeff[5] * p[1] * dir[0] - dir[2];//// b=(b1d0+b2d1+2b3p0d0+2b5p1d1+b4p0d1+b4p1d0-d2)
		real c = coeff[0] + coeff[1] * p[0] + coeff[2] * p[1] + coeff[3] * p[0] * p[0] + coeff[4] * p[1] * p[1] + coeff[5] * p[0] * p[1] - p[2];//// c=(b0+b1p0+b2p1+b3p0p0+b5p1p1+b4p0p1-p2)
		if (a < 0) { b = -b; c = -c; a = -a; }//// correct the sign
		//// solve at^2+bt+c=0
		bool is_valid = abs(a) > 1e-8;
		real delta = b * b - 4 * a * c; real t = 0;
		if (delta > 0) {
			real sqrt_delta = sqrt(delta);
			t = (-b - sqrt_delta) * 0.5 / a;
			//// choose the positive solution
			if (t < -v_r) t = (-b + sqrt_delta) * 0.5 / a;
			if (t < -v_r) is_valid = false;
			VectorD hit_p = p + t * dir;
			if (pow(hit_p(0), 2) + pow(hit_p(1), 2) > 0.25 * v_r * v_r) is_valid = false;
		}
		else is_valid = false;
		return is_valid ? t : -1e10;
	}
	return 0;
}

template<int d>
void PointSet<d>::Print_Statistics()
{
	using namespace AuxFunc;
	Seperation_Line(2);
	int avg_nb_num = 0;
	int min_nb_num = std::numeric_limits<int>::max();
	int max_nb_num = 0;
	int n = 0;
	for (int i = 0; i < points->Size(); i++) {
		if (points->I(i) == -1)continue;
		int i_nb_num = (int)tang_nbs[i].size();
		avg_nb_num += i_nb_num; n++;
		if (i_nb_num < min_nb_num) min_nb_num = i_nb_num;
		if (i_nb_num > max_nb_num) max_nb_num = i_nb_num;
	}
	if (n != 0)avg_nb_num /= n;
	std::cout << "[#Point] active: " << n << ", total: " << points->Size() << std::endl;
	std::cout << "[#Nb] avg: " << avg_nb_num << ", min: " << min_nb_num << ", max: " << max_nb_num << std::endl;
	std::cout << "[Nden] avg: " << avg_nden << ", init: " << init_avg_nden << std::endl;
	Seperation_Line(2);
}

template class PointSet<2>;
template class PointSet<3>;