#ifndef __Triangulation2D_h__
#define __Triangulation2D_h__
#include "Mesh.h"

namespace Triangulation{
bool Triangulation2D(const Array<Array<Vector2> >& loops,TriangleMesh<2>& mesh,const real max_edge_length);
};

#endif