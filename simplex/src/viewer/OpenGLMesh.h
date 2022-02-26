//////////////////////////////////////////////////////////////////////////
// Opengl mesh
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLObjectMesh_h__
#define __OpenGLObjectMesh_h__
#include <string>
#include "glm.hpp"
#include "gtc/matrix_transform.hpp"
#include "gtc/type_ptr.hpp"
#include "OpenGLColorRamp.h"
#include "Mesh.h"
#include "MeshFunc.h"
#include "GeometryPrimitives.h"
#include "OpenGLPrimitives.h"
#include "OpenGLObject.h"
#include "OpenGLTexture.h"
#include "OpenGLShaderLibrary.h"
#include "OpenGLShaderProgram.h"
#include "OpenGLFbos.h"
#include "OpenGLUbos.h"

const OpenGLColor default_mesh_color=OpenGLColor::Gray(.8f,1.f);

template<class T_MESH> class OpenGLMesh: public OpenGLObject
{public:typedef OpenGLObject Base;typedef T_MESH MESH_TYPE;using Base::Add_Shader_Program;
	T_MESH mesh;
	Field<Vector3,1> displacements;

	OpenGLMesh(){}

	virtual void Initialize()
	{
		Base::Initialize();
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vcolor"));
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vnormal_lt"));
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vtex"));
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("sd_depth"));
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("sd_lt"));
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vnormal_vtex_lt"));
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vnormal_vtex_lt_env"));
	}

	virtual void Update_Data_To_Render()
	{
		if(!Update_Data_To_Render_Pre())return;

		if(multi_color.size()>0){
			int n=(int)multi_color.size();
			int i=0;for(auto& p:mesh.Vertices()){
				OpenGL_Vertex4_And_Color4(Vertex(i),multi_color[i%n].rgba,opengl_vertices);i++;}}
		else{int i=0;for(auto& p:mesh.Vertices()){
			OpenGL_Vertex4_And_Color4(Vertex(i),color.rgba,opengl_vertices);i++;}}		////position, 4 floats; color, 4 floats		
		for(auto& e:mesh.elements)OpenGL_Vertex(e,opengl_elements);

		Set_OpenGL_Vertices();
		Set_OpenGL_Vertex_Attribute(0,4,8,0);	////position
		Set_OpenGL_Vertex_Attribute(1,4,8,4);	////color
		Set_OpenGL_Elements();

		Update_Data_To_Render_Post();
	}

	virtual void Refresh(const int frame)
	{
		bool is_binary_file=(File::File_Extension_Name(name)!="txt");
		std::string file_name=Object_File_Name(output_dir,frame,name);
		bool data_refreshed=false;
		if(is_binary_file){
			if(File::File_Exists(file_name)){
				mesh.elements.clear();
				File::Read_Binary_From_File(file_name,mesh);
				data_refreshed=true;
				if(verbose)std::cout<<"Read file "<<file_name<<std::endl;}}
		else{
			if(File::File_Exists(file_name)){
				mesh.elements.clear();
				File::Read_Text_From_File(file_name,mesh);
				data_refreshed=true;
				if(verbose)std::cout<<"Read file "<<file_name<<std::endl;}}

		if(use_vtx_displacement){
			file_name=Object_File_Name(output_dir,frame,displacement_name);
			if(File::File_Exists(file_name)){
				File::Read_Binary_From_File(file_name,displacements);
				data_refreshed=true;}
			else{
				std::cout<<"cannot read "<<displacement_name<<std::endl;
				displacements.Resize((int)mesh.Vertices().size(),Vector3::Zero());}}

		if(data_refreshed)Set_Data_Refreshed();
	}

	Vector3 Vertex(const int i) const 
	{return use_vtx_displacement?(mesh.Vertices()[i]+displacements.array[i]*scale):mesh.Vertices()[i];}
};

class OpenGLSegmentMesh : public OpenGLMesh<SegmentMesh<3> >
{public:typedef OpenGLMesh<SegmentMesh<3> > Base;
    OpenGLSegmentMesh(){color=default_mesh_color;name="segment_mesh";}

	virtual void Display() const
    {
		if(!visible||mesh.elements.empty())return;
		Update_Polygon_Mode();
		GLfloat old_line_width;glGetFloatv(GL_LINE_WIDTH,&old_line_width);
		{std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
		shader->Begin();
		glLineWidth(line_width);
		OpenGLUbos::Bind_Uniform_Block_To_Ubo(shader,"camera");
		glBindVertexArray(vao);
		glDrawElements(GL_LINES,ele_size,GL_UNSIGNED_INT,0);
		glLineWidth(old_line_width);
		shader->End();}
    }
};

class OpenGLTriangleMesh : public OpenGLMesh<TriangleMesh<3> >
{public:typedef OpenGLMesh<TriangleMesh<3> > Base;using Base::model;
	std::shared_ptr<OpenGLFbos::OpenGLFbo> fbo;
	glm::mat4 shadow_pv=glm::mat4(1.f);

    OpenGLTriangleMesh(){color=default_mesh_color;name="triangle_mesh";shading_mode=ShadingMode::Lighting;}

	virtual void Set_Shading_Mode(const ShadingMode& _mode)
	{
		shading_mode=_mode;
		switch(shading_mode){
		case ShadingMode::Shadow:{use_preprocess=true;use_depth_fbo=true;}break;}
	}

	virtual void Update_Data_To_Render()
	{
		if(!Update_Data_To_Render_Pre())return;

		switch(shading_mode){
		case ShadingMode::None:
		{use_vtx_color=true;use_vtx_normal=false;}break;
		case ShadingMode::Lighting:
		case ShadingMode::Shadow:
		case ShadingMode::EnvMapping:
		{use_vtx_color=false;use_vtx_normal=true;}break;
		case ShadingMode::TexOnly:
		{use_vtx_tex=true;use_vtx_color=false;use_vtx_normal=false;}break;
		case ShadingMode::TexLighting:
		{use_vtx_tex=true;use_vtx_color=false;use_vtx_normal=true;}break;}

		//if(mesh.normals==nullptr)use_vtx_normal=false;
		//else if((*mesh.normals).size()<mesh.Vertices().size())MeshFunc::Update_Normals(mesh);
		
		if(use_vtx_normal){
			if(mesh.normals==nullptr)mesh.normals.reset(new Array<Vector3>());
			if((*mesh.normals).size()<mesh.Vertices().size()||recomp_vtx_normal)MeshFunc::Update_Normals(mesh);}

		GLuint stride_size=4+(use_vtx_color?4:0)+(use_vtx_normal?4:0)+(use_vtx_tex?4:0);
		int i=0;for(auto& p:mesh.Vertices()){
			OpenGL_Vertex4(Vertex(i),opengl_vertices);	////position, 4 floats
			if(use_vtx_color){
				OpenGL_Color4(color.rgba,opengl_vertices);}				////color, 4 floats
			if(use_vtx_normal){
				OpenGL_Vertex4((*mesh.normals)[i],opengl_vertices);}	////normal, 4 floats
			if(use_vtx_tex){
				OpenGL_Vertex4((*mesh.textures)[i],opengl_vertices);}	////texture, 4 floats
			i++;}
		for(auto& e:mesh.elements)OpenGL_Vertex(e,opengl_elements);

		Set_OpenGL_Vertices();
		int idx=0;{Set_OpenGL_Vertex_Attribute(0,4,stride_size,0);idx++;}	////position
		if(use_vtx_color){Set_OpenGL_Vertex_Attribute(idx,4,stride_size,idx*4);idx++;}	////color
		if(use_vtx_normal){Set_OpenGL_Vertex_Attribute(idx,4,stride_size,idx*4);idx++;}	////normal
		if(use_vtx_tex){Set_OpenGL_Vertex_Attribute(idx,4,stride_size,idx*4);idx++;}////texture

		Set_OpenGL_Elements();

		Update_Data_To_Render_Post();
	}

	void Compute_Shadow_PV(glm::mat4& shadow_pv)	////change the viewport
	{
		using namespace OpenGLUbos;using namespace OpenGLFbos;
		Lights* lights=Get_Lights();if(lights==nullptr)return;
		Light* lt=lights->First_Shadow_Light();if(lt==nullptr)return;
		glm::vec3 light_pos=glm::vec3(lt->pos);

		fbo=Get_And_Bind_Fbo("depth",/*depth*/1);
		glViewport(0,0,fbo->width,fbo->height);

		glm::mat4 proj,view;
		float np=.001f,fp=10.f;
		float r=(float)fbo->height/(float)fbo->width;float w=5.f;float h=w*r;

		////TOFIX: fix bug for perspective
		//if(lt->Get_Type()==0)
			proj=glm::ortho(-w,w,-h,h,np,fp);
		//else{proj=glm::perspective(45.f*(float)deg_to_rad,1.f/r,np,fp);fbo->Set_Near_And_Far_Plane(np,fp);}

		view=glm::lookAt(light_pos,glm::vec3(0.f),glm::vec3(1.f,0.f,0.f));
		shadow_pv=proj*view;
	}

	virtual void Preprocess()
	{
		using namespace OpenGLUbos;using namespace OpenGLFbos;
		if(shading_mode!=ShadingMode::Shadow)return;

		GLint viewport[4];glGetIntegerv(GL_VIEWPORT,viewport);
		Compute_Shadow_PV(shadow_pv);

		std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[3];
		shader->Begin();
		//glEnable(GL_CULL_FACE);
		//glCullFace(GL_FRONT);
		Bind_Uniform_Block_To_Ubo(shader,"camera");
		shader->Set_Uniform_Matrix4f("shadow_pv",glm::value_ptr(shadow_pv));
		if(use_model){shader->Set_Uniform_Matrix4f("model",glm::value_ptr(model));}
		glBindVertexArray(vao);
		glDrawElements(GL_TRIANGLES,ele_size,GL_UNSIGNED_INT,0);
		//glDisable(GL_CULL_FACE);
		shader->End();
		Unbind_Fbo();

		glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
	}

	virtual void Display() const
    {
		using namespace OpenGLUbos;using namespace OpenGLFbos;using namespace OpenGLTextures;
		if(!visible||mesh.elements.empty())return;
		Update_Polygon_Mode();

		switch(shading_mode){
		case ShadingMode::None:{
			std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
			shader->Begin();
			Bind_Uniform_Block_To_Ubo(shader,"camera");
			glBindVertexArray(vao);
			glDrawElements(GL_TRIANGLES,ele_size,GL_UNSIGNED_INT,0);
			shader->End();
		}break;
		case ShadingMode::Lighting:{
			std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[1];
			shader->Begin();
			if(use_mat)shader->Set_Uniform_Mat(&mat);
			if(use_model){shader->Set_Uniform_Matrix4f("model",glm::value_ptr(model));}
			Bind_Uniform_Block_To_Ubo(shader,"camera");
			Bind_Uniform_Block_To_Ubo(shader,"lights");
			glBindVertexArray(vao);
			glDrawElements(GL_TRIANGLES,ele_size,GL_UNSIGNED_INT,0);
			shader->End();
		}break;
		case ShadingMode::Shadow:{
			std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[4];
			shader->Begin();
			if(use_mat)shader->Set_Uniform_Mat(&mat);
			Bind_Uniform_Block_To_Ubo(shader,"camera");
			Bind_Uniform_Block_To_Ubo(shader,"lights");
			if(fbo!=nullptr)fbo->Bind_As_Texture(0);
			shader->Set_Uniform("shadow_map",0);
			shader->Set_Uniform_Matrix4f("shadow_pv",glm::value_ptr(shadow_pv));
			if(use_model){shader->Set_Uniform_Matrix4f("model",glm::value_ptr(model));}
			glBindVertexArray(vao);
			glDrawElements(GL_TRIANGLES,ele_size,GL_UNSIGNED_INT,0);
			shader->End();
		}break;
		case ShadingMode::EnvMapping:{
			std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[6];
			shader->Begin();
			if(use_mat)shader->Set_Uniform_Mat(&mat);
			if(use_model){shader->Set_Uniform_Matrix4f("model",glm::value_ptr(model));}
			Bind_Uniform_Block_To_Ubo(shader,"camera");
			Bind_Uniform_Block_To_Ubo(shader,"lights");
			OpenGLTextures::Bind_Texture(env_name,TextureType::TxCube);
			shader->Set_Uniform("tex3d",0);
			glBindVertexArray(vao);
			glDrawElements(GL_TRIANGLES,ele_size,GL_UNSIGNED_INT,0);
			shader->End();
		}break;
		case ShadingMode::TexOnly:{
			std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[2];
			shader->Begin();
			Bind_Uniform_Block_To_Ubo(shader,"camera");
			if(use_vtx_tex)OpenGLTextures::Bind_Texture(tex_name,TextureType::Tx2d);
			shader->Set_Uniform("tex2d",0);
			glBindVertexArray(vao);
			glDrawElements(GL_TRIANGLES,ele_size,GL_UNSIGNED_INT,0);
			OpenGLTextures::Unbind_Texture();
			shader->End();
		}break;
		case ShadingMode::TexLighting:{
			std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[5];
			shader->Begin();
			if(use_mat)shader->Set_Uniform_Mat(&mat);
			if(use_model){shader->Set_Uniform_Matrix4f("model",glm::value_ptr(model));}
			Bind_Uniform_Block_To_Ubo(shader,"camera");
			Bind_Uniform_Block_To_Ubo(shader,"lights");
			if(use_vtx_tex)OpenGLTextures::Bind_Texture(tex_name,TextureType::Tx2d);
			shader->Set_Uniform("tex2d",0);
			glBindVertexArray(vao);
			glDrawElements(GL_TRIANGLES,ele_size,GL_UNSIGNED_INT,0);
			OpenGLTextures::Unbind_Texture();
			shader->End();
		}break;}
    }

	void Initialize_Cube_Mesh(const GLfloat scale=(real)1.,const Vector3& offset=Vector3::Zero(),bool inner_facing=false)
	{
		GLfloat vtx_pos[]={-1.f,1.f,-1.f,-1.f,-1.f,-1.f,1.f,-1.f,-1.f,1.f,-1.f,-1.f,1.f,1.f,-1.f,-1.f,1.f,-1.f,
			-1.f,-1.f,1.f,-1.f,-1.f,-1.f,-1.f,1.f,-1.f,-1.f,1.f,-1.f,-1.f,1.f,1.f,-1.f,-1.f,1.f,
			1.f,-1.f,-1.f,1.f,-1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,-1.f,1.f,-1.f,-1.f,
			-1.f,-1.f,1.f,-1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,1.f,-1.f,1.f,-1.f,-1.f,1.f,
			-1.f,1.f,-1.f,1.f,1.f,-1.f,1.f,1.f,1.f,1.f,1.f,1.f,-1.f,1.f,1.f,-1.f,1.f,-1.f,
			-1.f,-1.f,-1.f,-1.f,-1.f,1.f,1.f,-1.f,-1.f,1.f,-1.f,-1.f,-1.f,-1.f,1.f,1.f,-1.f,1.f};
		int vtx_num=36;mesh.Vertices().resize(vtx_num);for(int i=0;i<vtx_num;i++){mesh.Vertices()[i]=Vector3(vtx_pos[i*3],vtx_pos[i*3+1],vtx_pos[i*3+2])*scale+offset;}
		int ele_num=12;mesh.elements.resize(ele_num);for(int i=0;i<ele_num;i++)mesh.elements[i]=inner_facing?Vector3i(i*3,i*3+1,i*3+2):Vector3i(i*3,i*3+2,i*3+1);
	}

	void Initialize_Plane_Mesh(const real scale=(real)1.,const Vector3& offset=Vector3::Zero())
	{
		GLfloat vtx_pos[]={-1.f,0.f,-1.f, -1.f,0.f,1.f, 1.f,0.f,-1.f, -1.f,0.f,1.f, 1.f,0.f,1.f, 1.f,0.f,-1.f};
		GLfloat vtx_uv[]={0.f,1.f, 0.f,0.f, 1.f,1.f, 0.f,0.f, 1.f,0.f, 1.f,1.f};
		int vtx_num=6;mesh.Vertices().resize(vtx_num);for(int i=0;i<vtx_num;i++){mesh.Vertices()[i]=Vector3(vtx_pos[i*3],vtx_pos[i*3+1],vtx_pos[i*3+2])*scale+offset;}
		int ele_num=2;mesh.elements.resize(ele_num);for(int i=0;i<ele_num;i++)mesh.elements[i]=Vector3i(i*3,i*3+1,i*3+2);
		(*mesh.textures).resize(vtx_num);for(int i=0;i<vtx_num;i++)(*mesh.textures)[i]=Vector2(vtx_uv[i*2],vtx_uv[i*2+1]);
	}
};

class OpenGLTetrahedronMesh : public OpenGLMesh<TetrahedronMesh<3> >
{public:typedef OpenGLMesh<TetrahedronMesh<3> > Base;
    OpenGLTetrahedronMesh(){color=default_mesh_color;name="tetrahedron_mesh";polygon_mode=PolygonMode::Wireframe;}

	virtual void Update_Data_To_Render()
	{
		if(!Update_Data_To_Render_Pre())return;

		int i=0;for(auto& p:mesh.Vertices()){
			OpenGL_Vertex4_And_Color4(Vertex(i),color.rgba,opengl_vertices);i++;}	////position, 4 floats; color, 4 floats

		switch(polygon_mode){
		case PolygonMode::Wireframe:{
			for(auto& e:mesh.elements){
				OpenGL_Vertex(e[0],e[1],e[2],opengl_elements);
				OpenGL_Vertex(e[0],e[2],e[3],opengl_elements);
				OpenGL_Vertex(e[0],e[3],e[1],opengl_elements);
				OpenGL_Vertex(e[1],e[2],e[3],opengl_elements);}
		}break;
		case PolygonMode::Fill:case PolygonMode::SurfOnly:{
			Array<Vector3i> surf_elements;MeshFunc::Volumetric_Mesh_Surface(mesh.elements,surf_elements);
			for(auto& e:surf_elements)OpenGL_Vertex(e,opengl_elements);
		}break;}

		Set_OpenGL_Vertices();
		Set_OpenGL_Vertex_Attribute(0,4,8,0);	////position
		Set_OpenGL_Vertex_Attribute(1,4,8,4);	////color
		Set_OpenGL_Elements();

		Update_Data_To_Render_Post();
	}

	virtual void Display() const
    {
    	using namespace OpenGLUbos;
		if(!visible||mesh.elements.empty())return;
		Update_Polygon_Mode();

		{std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
		shader->Begin();
		Bind_Uniform_Block_To_Ubo(shader,"camera");
		glBindVertexArray(vao);
		glDrawElements(GL_TRIANGLES,ele_size,GL_UNSIGNED_INT,0);
		shader->End();}
    }

protected:
	Array<Vector3i> surface_elements;
};

class OpenGLQuadMesh : public OpenGLMesh<QuadMesh<3> >
{public:typedef OpenGLMesh<QuadMesh<3> > Base;
    OpenGLQuadMesh(){color=default_mesh_color;name="quad_mesh";polygon_mode=PolygonMode::Wireframe;}

	virtual void Update_Data_To_Render()
	{
		if(!Update_Data_To_Render_Pre())return;

		int i=0;for(auto& p:mesh.Vertices()){
			OpenGL_Vertex4_And_Color4(Vertex(i),color.rgba,opengl_vertices);i++;}	////position, 4 floats; color, 4 floats

		switch(polygon_mode){
		case PolygonMode::Wireframe:{
			for(auto e:mesh.elements){e={e[0],e[2],e[3],e[1]};	////correct the default order in an quad element, making the four vertices follow counterclockwise order
				OpenGL_Vertex(e,opengl_elements);opengl_elements.push_back(UINT_MAX);}	////quad, 5 uint
		}break;
		case PolygonMode::Fill:case PolygonMode::SurfOnly:{std::cerr<<"Error: [OpenGLQuadMesh] draw mode not implemented"<<std::endl;}break;}	////TODO

		Set_OpenGL_Vertices();
		Set_OpenGL_Vertex_Attribute(0,4,8,0);	////position
		Set_OpenGL_Vertex_Attribute(1,4,8,4);	////color
		Set_OpenGL_Elements();

		Update_Data_To_Render_Post();
	}

	virtual void Display() const
    {
		using namespace OpenGLUbos;
		if(!visible||mesh.elements.empty())return;
		Update_Polygon_Mode();

		{std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
		shader->Begin();
		Bind_Uniform_Block_To_Ubo(shader,"camera");
		glBindVertexArray(vao);
		glEnable(GL_PRIMITIVE_RESTART);
		glPrimitiveRestartIndex((GLuint)UINT_MAX);
		glDrawElements(GL_LINE_LOOP,ele_size,GL_UNSIGNED_INT,0);
		glDisable(GL_PRIMITIVE_RESTART);
		shader->End();}
    }
};

class OpenGLSkybox : public OpenGLTriangleMesh
{public:typedef OpenGLTriangleMesh Base;
	GLfloat scale;
	std::string tex_name;

    OpenGLSkybox(){color=default_mesh_color;name="skybox";env_name="skybox";scale=(real)10.;}

	virtual void Initialize()
	{
		OpenGLObject::Initialize();
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("skybox"));
		Initialize_Cube_Mesh(scale);
	}

	virtual void Update_Data_To_Render()
	{
		if(!Update_Data_To_Render_Pre())return;

		for(auto& p:mesh.Vertices()){
			OpenGL_Vertex4(p,opengl_vertices);}					////position, 4 floats
		for(auto& e:mesh.elements)OpenGL_Vertex(e,opengl_elements);

		Set_OpenGL_Vertices();
		Set_OpenGL_Vertex_Attribute(0,4,4,0);	////position
		Set_OpenGL_Elements();
		Update_Data_To_Render_Post();
	}

	virtual void Display() const
    {
		using namespace OpenGLUbos;using namespace OpenGLTextures;
		if(!visible)return;
		Update_Polygon_Mode();

		std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
		shader->Begin();
		Bind_Uniform_Block_To_Ubo(shader,"camera");
		OpenGLTextures::Bind_Texture(env_name,TextureType::TxCube);
		shader->Set_Uniform("tex3d",0);
		glBindVertexArray(vao);
		glDrawElements(GL_TRIANGLES,ele_size,GL_UNSIGNED_INT,0);
		OpenGLTextures::Unbind_Texture();
		shader->End();		
    }
};
#endif