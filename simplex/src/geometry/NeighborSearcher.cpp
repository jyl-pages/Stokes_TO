//////////////////////////////////////////////////////////////////////////
// Implementations for NeighborSearcher
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "NeighborSearcher.h"

template<int d>
void NeighborSearcher<d>::Update_Points(Array<VectorD>& points)
{
	search_results.clear();
	this->Build_Data(points);
}

template<int d>
void NeighborSearcher<d>::Update_Points(Array<VectorD>& points, FilterFunc& filter_func)
{
	Array<VectorD> temp_array;
	temp_array.clear(); temp_array.reserve(points.size());
	for (size_t i = 0; i < points.size(); i++) {
		if (filter_func((int)i)) temp_array.push_back(points[i]);
	}
	this->Update_Points(temp_array);
}

template<int d>
size_t NeighborSearcher<d>::Find_Neighbors(const VectorD& pos, const real& radius, FilterFunc& filter_func, Array<int>& results, bool append) const
{
	if (!append) results.clear();
	Array<int> temp_results;
	this->Find_Neighbors(pos, radius, temp_results, false);//append=false, temp_results is cleared
	size_t num = 0;
	for (size_t i = 0; i < temp_results.size(); i++) {
		if (filter_func(temp_results[i])) {
			num++;
			results.push_back(temp_results[i]);
		}
	}
	return num;
}

template<int d>
Array<int> NeighborSearcher<d>::Find_Neighbors(const VectorD& pos, const real& radius) const
{
	Array<int> temp_results;
	this->Find_Neighbors(pos, radius, temp_results, false);//append=false, temp_results is cleared
	return temp_results;
}

template<int d>
Array<int> NeighborSearcher<d>::Find_Neighbors(const VectorD& pos, const real& radius, FilterFunc& filter_func) const
{
	Array<int> temp_results;
	this->Find_Neighbors(pos, radius, filter_func, temp_results, false);//append=false, temp_results are cleared
	return temp_results;
}

template<int d>
ArraySlice<int> NeighborSearcher<d>::Record_Neighbors(const VectorD& pos, const real& radius)
{
	size_t beg = search_results.size();
	this->Find_Neighbors(pos, radius, search_results, true);//append=true
	size_t end = search_results.size();
	return ArraySlice<int>(beg, end, search_results);
}

template<int d>
ArraySlice<int> NeighborSearcher<d>::Record_Neighbors(const VectorD& pos, const real& radius, FilterFunc& filter_func)
{
	size_t beg = search_results.size();
	this->Find_Neighbors(pos, radius, filter_func, search_results, true);//append=true
	size_t end = search_results.size();
	return ArraySlice<int>(beg, end, search_results);
}


template<int d>
void PointSetAdapter<d>::Initialize(Array<VectorD>& arr) {
	points_ptr = std::make_shared<Array<VectorD> >(arr);
}

template<int d>
size_t PointSetAdapter<d>::kdtree_get_point_count()const
{
	if (!points_ptr) return 0;//if it's null
	return points_ptr->size();
}

template<int d>
real PointSetAdapter<d>::kdtree_get_pt(const size_t idx, const size_t dim) const
{
	if (!points_ptr) return 0;
	return (*points_ptr)[idx][dim];
}

template<int d>
void NeighborKDTree<d>::Build_Data(Array<VectorD>& arr)
{
	points.Initialize(arr);
	index.buildIndex();
}

template<int d>
size_t NeighborKDTree<d>::Find_Neighbors(const VectorD& pos, const real& radius, Array<int>& results, bool append) const
{
	nanoflann::SearchParams params; params.sorted = false;
	std::vector<std::pair<size_t, real> > ret_matches;//it will be cleared in radiusSearch()
	size_t nMatches = index.radiusSearch(pos.data(), radius * radius, ret_matches, params);//L2 is squared
	if (!append) results.clear();
	for (size_t i = 0; i < ret_matches.size(); i++) {
		results.push_back((int)ret_matches[i].first);
	}
	return nMatches;
}

template<int d>
int NeighborKDTree<d>::Find_Nearest_Nb(const VectorD& pos) const
{
	std::array<size_t, 1> ret_index;
	std::array<real, 1> out_dist_sqr;
	size_t num_results = index.knnSearch(pos.data(), 1, &ret_index[0], &out_dist_sqr[0]);
	if (num_results == 0) {
		std::cerr << "[Error]NeighborKDTree::Find_Nearest_Nb fails for " << pos.transpose() << std::endl;
		exit(1);
		return -1;
	}
	else {
		return (int)ret_index[0];
	}
}

template<int d>
int NeighborKDTree<d>::Find_K_Nearest_Nb(const VectorD& pos, int k, Array<int>& results)const
{
	size_t* ret_index = new size_t[k];
	real* out_dist_sqr = new real[k];
	size_t num_results = index.knnSearch(pos.data(), k, &ret_index[0], &out_dist_sqr[0]);
	results.resize(num_results);
	for (size_t i = 0; i < num_results; i++) {
		size_t ith = ret_index[i];
		results[i] = (int)ith;
	}
	delete[] ret_index;
	delete[] out_dist_sqr;
	return (int)num_results;
}

template class NeighborSearcher<2>;
template class NeighborSearcher<3>;

template class PointSetAdapter<2>;
template class PointSetAdapter<3>;

template class NeighborKDTree<2>;
template class NeighborKDTree<3>;
