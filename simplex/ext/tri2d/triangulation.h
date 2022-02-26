#ifndef __triangulation_h__
#define __triangulation_h__

#include <iostream>

////////////////////////////////
int delaunay_triangulation(
						   // input
						   int nloop,	// number of loops (at least one = outmost loop)  e.g 2 (dounat)
						   int* nxys,	// numbers of vertices in each loop (first one is the outmost loop)  e.g. {4, 3} (triangle inside square)
						   double* xys,	// xyz values  e.g. {0,0,  0,1, 1,1, 1,0,      0.1,0.1, 0,5,0.9, 0.9,0.1} 
						   double max_edge_length, // maximum edge length (if minus then no adding) e.g. {-1}

						   // output
						   int* ntri,		// number of triangles e.g. 2
						   int** atri,		// list of triangles (each triangle is a list of three indices) e.g. {0,1,2, 0,2,3}
						   int* nxys_new,	// number of newly added vertices e.g. 1
						   double** xys_new);	// coordinates of newly added vertices e.g. {0,0.5,0}

void test_dtri();

#endif