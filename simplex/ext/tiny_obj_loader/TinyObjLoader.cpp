//#####################################################################
// Tiny obj loader
// Bo Zhu
//#####################################################################
#include <iostream>
#include <memory>
#include "Common.h"
#include "Mesh.h"
#include "tiny_obj_loader.h"
#include "TinyObjLoader.h"

namespace Obj{

template<class T_MESH> void Read_From_Obj_File(const std::string& file_name,Array<std::shared_ptr<T_MESH> >& meshes)
{
    Array<tinyobj::shape_t> shapes;Array<tinyobj::material_t> materials;std::string base;size_t l;
    if ((l=file_name.find_last_of('/')) != std::string::npos)base=file_name.substr(0, l+1);
    else if((l=file_name.find_last_of('\\'))!=std::string::npos)base=file_name.substr(0, l+1);
    std::string err=tinyobj::LoadObj(shapes,materials,file_name.c_str(),base.c_str());
    if(err!="")std::cerr<<"Read obj error: "<<err<<std::endl;
    else std::cout<<"Read obj file: "<<file_name<<", #shapes="<<shapes.size()<<", #materials="<<materials.size()<<std::endl;

    meshes.resize((int)shapes.size());
    for(auto i=0;i<shapes.size();i++){meshes[i]=std::make_shared<T_MESH>();tinyobj::mesh_t& mesh=shapes[i].mesh;
        for(auto j=0;j<mesh.positions.size()/3;j++){
            (*meshes[i]->vertices).push_back(Vector3((real)mesh.positions[j*3],(real)mesh.positions[j*3+1],(real)mesh.positions[j*3+2]));}
        for(auto j=0;j<mesh.indices.size()/3;j++){
            meshes[i]->elements.push_back(Vector3i(mesh.indices[j*3],mesh.indices[j*3+1],mesh.indices[j*3+2]));}
        if (!meshes[i]->normals) meshes[i]->normals = std::make_shared< Array<Vector3> >();
        for (auto j=0; j<mesh.normals.size()/3; j++) {
            (*meshes[i]->normals).push_back(Vector3((real)mesh.normals[j*3],(real)mesh.normals[j*3+1],(real)mesh.normals[j*3+2]));}
    }
}

template void Read_From_Obj_File<TriangleMesh<3> >(const std::string&,Array<std::shared_ptr<TriangleMesh<3> > >&);

};