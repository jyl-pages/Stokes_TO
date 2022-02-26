#include "Triangulation2D.h"
#include "triangulation.h"

namespace Triangulation
{
bool Triangulation2D(const Array<Array<Vector2> >& loops,TriangleMesh<2>& mesh,const real max_edge_length)
{
#ifdef USE_TRI2D
	int loop_n=(int)loops.size();int vtx_n=0;
	Array<int> vtx_n_in_loop(loop_n,0);
	for(int i=0;i<loop_n;i++){vtx_n_in_loop[i]=(int)loops[i].size();vtx_n+=vtx_n_in_loop[i];}

	Array<real> pos(vtx_n*2,(real)0);int p=0;
	for(int i=0;i<loop_n;i++){
		for(int j=0;j<vtx_n_in_loop[i];j++){
			pos[p*2]=loops[i][j][0];
			pos[p*2+1]=loops[i][j][1];p++;}}

	int tri_n=0;
	int* tri=nullptr;
	int vtx_n_new=0;
	real* pos_new=nullptr;

	bool result=delaunay_triangulation(loop_n,&vtx_n_in_loop[0],&pos[0],max_edge_length,&tri_n,&tri,&vtx_n_new,&pos_new);
	if(!result)return false;

	mesh.Vertices().resize(vtx_n+vtx_n_new);
	for(int i=0;i<vtx_n;i++){mesh.Vertices()[i]=Vector2(pos[i*2],pos[i*2+1]);}
	for(int i=0;i<vtx_n_new;i++){mesh.Vertices()[vtx_n+i]=Vector2(pos_new[i*2],pos_new[i*2+1]);}
	mesh.Elements().resize(tri_n);
	for(int i=0;i<tri_n;i++){mesh.Elements()[i]=Vector3i(tri[i*3],tri[i*3+1],tri[i*3+2]);}
	std::cout<<"Delaunay triangulation mesh: "<<mesh.Vertices().size()<<", "<<mesh.Elements().size()<<std::endl;

	return true;
#else
	std::cerr<<"Error: [Triangulation2D] requires USE_TRI2D"<<std::endl;
	return false;
#endif
}
};