//#####################################################################
// TetGen Func
// Jingyuan Zhu
//#####################################################################
#ifndef __TetGenFunc_h__
#define __TetGenFunc_h__
#include <string>
#include <memory>
#include <vector>
#include "Mesh.h"

namespace TetGen
{
    void Generate(const TriangleMesh<3>& surface_mesh,TetrahedronMesh<3>& volume_mesh,const real dx);
};

#endif