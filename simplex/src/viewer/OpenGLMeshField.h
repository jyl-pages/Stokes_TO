//////////////////////////////////////////////////////////////////////////
// Opengl mesh field
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLMeshField_h__
#define __OpenGLMeshField_h__
#include "OpenGLField.h"

template<class T_VAL,class T_MESH,int data_size>
class OpenGLMeshField : public OpenGLField<Field<T_VAL,1>,data_size>
{typedef OpenGLField<Field<T_VAL,1>,data_size> Base;
public:
	using Base::color_mapper;using Base::field_val;using Base::data_idx;using Base::opengl_vertices;using Base::opengl_elements;using Base::name;
	using Base::visible;using Base::shader_programs;using Base::vao;using Base::ele_size;using Base::data;using Base::output_dir;using Base::alpha;
	using Base::field_draw_type;using Base::data_refreshed;using Base::polygon_mode;
	using Base::Update_Data_To_Render_Pre;using Base::Update_Data_To_Render_Post;using Base::Scalar;using Base::Set_OpenGL_Vertices;using Base::Set_OpenGL_Vertex_Attribute;using Base::Set_Data_Refreshed;
	using Base::Set_OpenGL_Elements;using Base::Use_Alpha_Blend;using Base::Disable_Alpha_Blend;using Base::Enable_Alpha_Blend;using Base::Update_Polygon_Mode;
	using Base::Read_Mesh_Field;using Base::Update_Data_To_Render_Mesh_Lines;

	const OpenGLMesh<T_MESH>* opengl_mesh;

	OpenGLMeshField(const OpenGLMesh<T_MESH>* _opengl_mesh=nullptr):Base(),opengl_mesh(_opengl_mesh){}
	void Set_Data_Pointers(const OpenGLMesh<T_MESH>* _opengl_mesh){opengl_mesh=_opengl_mesh;}

	virtual void Update_Data_To_Render()
	{
		if(!Update_Data_To_Render_Pre())return;

		if(opengl_mesh==nullptr){return;}
		auto& mesh=opengl_mesh->mesh;

		for(size_t i=0;i<mesh.Vertices().size();i++){
			const Vector1i node=Vec1i((int)i);
			const Vector3& pos=opengl_mesh->Vertex((int)i);
			OpenGLColor color=color_mapper->Color(Scalar((*field_val[data_idx])(node)),(int)data[data_idx].color_type);
			color.rgba[3]=alpha;
			OpenGL_Vertex4_And_Color4(pos,color.rgba,opengl_vertices);}		////position, 4 floats; color, 4 floats
		
		if(mesh.elements.size()>0){
			const int ele_dim=(int)mesh.elements[0].size();
			switch(ele_dim){
			case 3:{	////triangles
				for(auto& e:mesh.elements)OpenGL_Vertex(e,opengl_elements);
			}break;
			case 4:{	////tetrahedra
				for(auto& e:mesh.elements){
					OpenGL_Vertex(e[0],e[1],e[2],opengl_elements);
					OpenGL_Vertex(e[0],e[2],e[3],opengl_elements);
					OpenGL_Vertex(e[0],e[3],e[1],opengl_elements);
					OpenGL_Vertex(e[1],e[2],e[3],opengl_elements);}
			}break;}}

		Set_OpenGL_Vertices();
		Set_OpenGL_Vertex_Attribute(0,4,8,0);	////position
		Set_OpenGL_Vertex_Attribute(1,4,8,4);	////color
		Set_OpenGL_Elements();
		
		Update_Data_To_Render_Post();
	}

    virtual void Display() const
    {
		if(!visible)return;
		Update_Polygon_Mode();
		{std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
		shader->Begin();
		glBindVertexArray(vao);
		if(Use_Alpha_Blend())Enable_Alpha_Blend();
		glDrawElements(GL_TRIANGLES,ele_size,GL_UNSIGNED_INT,0);
		if(Use_Alpha_Blend())Disable_Alpha_Blend();
		shader->End();}
    }
	
	virtual void Refresh(const int frame)
	{
		Set_Data_Refreshed(false);	////TOFIX: why set false first?
		if(opengl_mesh==nullptr)return;
		for(size_type i=0;i<data.size();i++){std::string file_name=output_dir+"/"+std::to_string(frame)+"/"+data[i].name;
			if(File::File_Exists(file_name)){
				Read_Mesh_Field(file_name,(int)opengl_mesh->mesh.elements.size(),(int)opengl_mesh->mesh.vertices->size(),
					opengl_mesh->mesh.elements,field_val[i],data[i].store_type);Set_Data_Refreshed();}}
	}
};
#endif