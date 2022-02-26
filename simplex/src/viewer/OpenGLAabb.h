//////////////////////////////////////////////////////////////////////////
// Opengl AABB
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLAabb_h__
#define __OpenGLAabb_h__
#include "File.h"
#include "GeometryPrimitives.h"
#include "OpenGLObject.h"
#include "OpenGLShaderLibrary.h"
#include "OpenGLPrimitives.h"
#include "OpenGLUbos.h"
#include "OpenGLTexture.h"
#include "OpenGLFbos.h"
#include "OpenGLColorMapper.h"

using namespace OpenGLUbos;
using namespace OpenGLTextures;
using namespace OpenGLFbos;

class OpenGLAabb : public OpenGLObject
{public:typedef OpenGLObject Base;
	Box<3> aabb;
	
	OpenGLAabb():aabb(Vector3::Zero(),Vector3::Ones()){color=OpenGLColor::Gray(.2f,1.f);name="aabb";}

	virtual void Initialize()
	{
		OpenGLObject::Initialize();
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("aabb_ray_depth"));
	}

	virtual void Update_Data_To_Render()
	{
		if(!Update_Data_To_Render_Pre())return;
		Grid<3> cube_grid;cube_grid.Initialize(Vector3i::Ones());
		iterate_node_d(iter,cube_grid,3){const Vector3i& node=iter.Coord();
			Vector3 pos=aabb.min_corner+aabb.Edge_Lengths().cwiseProduct(node.cast<real>());
			OpenGL_Vertex4(pos,opengl_vertices);}	////position, 4 floats
		ArrayF2P<int,3> nodes;for(int i=0;i<nodes.size();i++)nodes[i]=i;
		OpenGL_Cube(nodes,opengl_elements);/*cube, 18 uint*/
		
		Set_OpenGL_Vertices();
		Set_OpenGL_Vertex_Attribute(0,4,4,0);	////position
		Set_OpenGL_Elements();	
		Update_Data_To_Render_Post();
	}

    virtual void Display() const
    {
		if(!visible)return;
		Update_Polygon_Mode();
		{std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
		shader->Begin();
		glEnable(GL_CULL_FACE);
		glCullFace(GL_FRONT);
		glBindVertexArray(vao);
		glEnable(GL_PRIMITIVE_RESTART);
		glPrimitiveRestartIndex((GLuint)UINT_MAX);
		glDrawElements(GL_TRIANGLE_STRIP,ele_size,GL_UNSIGNED_INT,0);
		glDisable(GL_PRIMITIVE_RESTART);
		glDisable(GL_CULL_FACE);
		shader->End();}
    }
};

////The template parameter is only for the data type. The texture data type is always GLfloat.
template<class T> class OpenGLVolume : public OpenGLAabb
{public:typedef OpenGLAabb Base;
	OpenGLAabb opengl_aabb;
	Vector2f screen_size;
	std::shared_ptr<OpenGLFbo> fbo;
	std::shared_ptr<OpenGLTextureVolume<> > tex_vol;
	std::shared_ptr<OpenGLTexture1D<ushort> > tex_trsf;
	GLfloat dx;
	Grid<3> grid;
	OpenGLColorTransfer c_trsf;

	OpenGLVolume():fbo(nullptr),tex_vol(nullptr),tex_trsf(nullptr),dx(.05f){color=OpenGLColor::Gray(.2f,1.f);name="volume";use_preprocess=true;}

	void Set_Data_Pointers(const Grid<3>* _grid,const Field<GLfloat,3>* field,const Array<ushort>* transfer=nullptr)
	{
		grid=*_grid;
		tex_vol->Update_Data_To_Render(field);
		tex_trsf->Update_Data_To_Render(transfer);
		opengl_aabb.aabb=Box<3>(grid.domain_min,grid.domain_max);
		aabb=opengl_aabb.aabb;
		dx=.75f*(GLfloat)grid.dx;
		Set_Data_Refreshed();
	}

	virtual void Initialize()
	{
		OpenGLObject::Initialize();
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vol_ray_casting"));

		tex_vol.reset(dynamic_cast<OpenGLTextureVolume<>*>(Get_Texture("volume",TextureType::Tx3d).get()));
		tex_trsf.reset(dynamic_cast<OpenGLTexture1D<ushort>*>(Get_Texture("1d",TextureType::Tx1d).get()));

		////default color transfer function
		////levelset
		//Array<ArrayF<GLfloat,5> > c_array={{.0f,0.f,0.f,0.f,0.f},{0.49f,0.f,0.f,0.f,0.f},{.5f,1.f,0.f,0.f,.01f},{.51f,0.f,0.f,0.f,0.f},{1.f,0.f,0.f,0.f,0.f}};
		////density field
		Array<ArrayF<GLfloat,5> > c_array={{0.f,0.f,0.f,0.f,0.f},{0.02f,0.f,0.f,0.f,0.f},{1.f,1.f,1.f,1.f,.9f}};
		c_trsf.Initialize(c_array);
	}

	virtual void Preprocess()
	{
		fbo=Get_And_Bind_Fbo("screen",/*position*/2);if(fbo==nullptr)return;

		fbo->Resize_To_Window();
		screen_size=Vector2f((GLfloat)fbo->width-.5f,(GLfloat)fbo->height-.5f);
		
		fbo->Clear();
		fbo->Bind();
		opengl_aabb.Display();
		fbo->Write_To_File("test");
		fbo->Unbind();
	}

	virtual void Update_Data_To_Render()
	{
		opengl_aabb.Update_Data_To_Render();
		Base::Update_Data_To_Render();
	}

	virtual void Set_Data_Refreshed(const bool _refreshed=true)
	{
		data_refreshed=_refreshed;
		opengl_aabb.data_refreshed=_refreshed;
	}

    virtual void Display() const
    {
		if(!visible)return;
		Vector3f org=opengl_aabb.aabb.min_corner.cast<float>();
		Vector3 length=opengl_aabb.aabb.Edge_Lengths();
		Vector3f one_over_length=length.cast<float>().cwiseInverse();

		Update_Polygon_Mode();
		{std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
		shader->Begin();
		Enable_Alpha_Blend();
		shader->Set_Uniform("screen_size",screen_size);
		shader->Set_Uniform("dx",dx);
		shader->Set_Uniform("org",org);
		shader->Set_Uniform("one_over_length",one_over_length);
		if(fbo!=nullptr)fbo->Bind_As_Texture(0);
		shader->Set_Uniform("tex2d",0);
		tex_vol->Bind(1);
		shader->Set_Uniform("tex3d",1);
		tex_trsf->Bind(2);
		shader->Set_Uniform("tex1d",2);
		glBindVertexArray(vao);
		glEnable(GL_CULL_FACE);
		glCullFace(GL_BACK);
		glEnable(GL_PRIMITIVE_RESTART);
		glPrimitiveRestartIndex((GLuint)UINT_MAX);
		glDrawElements(GL_TRIANGLE_STRIP,ele_size,GL_UNSIGNED_INT,0);
		glDisable(GL_PRIMITIVE_RESTART);
		glDisable(GL_CULL_FACE);
		Disable_Alpha_Blend();
		shader->End();}
    }

	virtual void Refresh(const int frame)
	{
		if(!initialized)Initialize();

		if(frame==0){std::string file_name=Object_File_Name(output_dir,frame,name);
			if(File::File_Exists(file_name)){File::Read_Binary_From_File(file_name,grid);}}
		
		if(data.size()>0){
			std::string file_name=Object_File_Name(output_dir,frame,data[0].name);
			if(File::File_Exists(file_name)){
				std::shared_ptr<Field<T,3> > read_field=std::make_shared<Field<T,3> >();
				File::Read_Binary_From_File(file_name,*read_field);
				Field<GLfloat,3> field;field.Resize(read_field->counts);

				iterate_cell_d(iter,grid,3){const Vector3i& cell=iter.Coord();
					int idx=field.counts[0]*field.counts[1]*cell[2]+field.counts[0]*cell[1]+cell[0];field.array[idx]=(GLfloat)(*read_field)(cell);}
				Set_Data_Pointers(&grid,&field,&c_trsf.transfer);
				Set_Data_Refreshed();}}
	}
};

template<class T,int d> class OpenGLVectorVolume : public OpenGLAabb
{public:typedef OpenGLAabb Base;
	OpenGLAabb opengl_aabb;
	Vector2f screen_size;
	std::shared_ptr<OpenGLFbo> fbo;
	std::shared_ptr<OpenGLTextureVectorVolume<T,d> > tex_vec_vol;
	ArrayF<std::shared_ptr<OpenGLTexture1D<T> >,d> tex_trsf;
	GLfloat dx;
	Grid<3> grid;
	glm::vec4 pl_pos;
	ArrayF<OpenGLColorTransfer,d> c_trsf;

	OpenGLVectorVolume():fbo(nullptr),tex_vec_vol(nullptr),dx(.05f){color=OpenGLColor::Gray(.2f,1.f);name="vector_volume";use_preprocess=true;}
	
	void Set_Data_Pointers(const Grid<3>* _grid,const Field<Vector<T,d>,3>* field,const ArrayF<Array<T>,d>* transfer=nullptr)
	{
		grid=*_grid;
		tex_vec_vol->Update_Data_To_Render(field);
		for(int i=0;i<d;i++)tex_trsf[i]->Update_Data_To_Render(&(*transfer)[i]);
		opengl_aabb.aabb=Box<3>(grid.domain_min,grid.domain_max);
		aabb=opengl_aabb.aabb;
	}

	void Set_Ray_Casting_Dx(const GLfloat _dx){dx=_dx;}

	virtual void Initialize()
	{
		OpenGLObject::Initialize();
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vec_vol_ray_casting"));

		std::string tex_name="vec"+std::to_string(d)+"_vol";
		TextureType type=(TextureType)((int)TextureType::Tx3d2+d-2);
		tex_vec_vol.reset(dynamic_cast<OpenGLTextureVectorVolume<T,d>*>(Get_Texture(tex_name,type).get()));
		for(int i=0;i<d;i++){std::string tex_name="1d";if(i>0)tex_name+=std::to_string(i+1);
			tex_trsf[i].reset(dynamic_cast<OpenGLTexture1D<T>*>(Get_Texture(tex_name,TextureType::Tx1d).get()));}

		GLfloat alpha=.05f;
		const ArrayF<GLfloat,4> alpha_coef={1.f,1.f,1.f,1.f};
		//const ArrayF<Vector3f,4> color={Vector3f(1.f,0.f,0.f),Vector3f(0.f,1.f,0.f),Vector3f(0.f,0.f,1.f),Vector3f(0.f,0.f,0.f)};
		const ArrayF<Vector3f,4> color={Vector3f(0.f,1.f,1.f),Vector3f(1.f,1.f,0.f),Vector3f(1.f,0.f,0.f),Vector3f(0.f,0.f,1.f)};

		for(int i=0;i<d;i++){
			switch(i){
			case 0:{Array<ArrayF<GLfloat,5> > c_array={{.0f,0.f,0.f,0.f,0.f},{.0001f,0.f,0.f,0.f,0.f},{.7f,color[i][0],color[i][1],color[i][2],alpha*alpha_coef[i]},{1.f,color[i][0],color[i][1],color[i][2],1.f*alpha_coef[i]}};
				c_trsf[i].Initialize(c_array);}break;
			case 1:{Array<ArrayF<GLfloat,5> > c_array={{.0f,0.f,0.f,0.f,0.f},{.0001f,0.f,0.f,0.f,0.f},{.7f,color[i][0],color[i][1],color[i][2],alpha*alpha_coef[i]},{1.f,color[i][0],color[i][1],color[i][2],1.f*alpha_coef[i]}};
				c_trsf[i].Initialize(c_array);}break;
			case 2:{Array<ArrayF<GLfloat,5> > c_array={{.0f,0.f,0.f,0.f,0.f},{.0001f,0.f,0.f,0.f,0.f},{.7f,color[i][0],color[i][1],color[i][2],alpha*alpha_coef[i]},{1.f,color[i][0],color[i][1],color[i][2],1.f*alpha_coef[i]}};
				c_trsf[i].Initialize(c_array);}break;
			case 3:{Array<ArrayF<GLfloat,5> > c_array={{.0f,0.f,0.f,0.f,0.f},{.0001f,0.f,0.f,0.f,0.f},{.7f,color[i][0],color[i][1],color[i][2],alpha*alpha_coef[i]},{1.f,color[i][0],color[i][1],color[i][2],1.f*alpha_coef[i]}};
				c_trsf[i].Initialize(c_array);}break;}}
	}

	virtual void Preprocess()
	{
		Camera* camera=Get_Camera();
		if(camera!=nullptr){pl_pos=camera->position;}

		fbo=Get_And_Bind_Fbo("screen",/*position*/2);if(fbo==nullptr)return;
		fbo->Resize_To_Window();
		screen_size=Vector2f((GLfloat)fbo->width-.5f,(GLfloat)fbo->height-.5f);
		
		fbo->Clear();
		fbo->Bind();
		opengl_aabb.Display();
		fbo->Unbind();
	}

	virtual void Update_Data_To_Render()
	{
		opengl_aabb.Update_Data_To_Render();
		Base::Update_Data_To_Render();
	}

    virtual void Display() const
    {
		if(!visible)return;
		Vector3f org=opengl_aabb.aabb.min_corner.cast<float>();
		Vector3 length=opengl_aabb.aabb.Edge_Lengths();
		Vector3f one_over_length=length.cast<float>().cwiseInverse();

		Update_Polygon_Mode();
		{std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
		shader->Begin();
		Enable_Alpha_Blend();
		shader->Set_Uniform("screen_size",screen_size);
		shader->Set_Uniform("dx",dx);
		shader->Set_Uniform("tex_dim",d);
		shader->Set_Uniform("org",org);
		shader->Set_Uniform("one_over_length",one_over_length);
		shader->Set_Uniform("pl_pos",pl_pos);
		if(fbo!=nullptr)fbo->Bind_As_Texture(0);
		shader->Set_Uniform("tex2d",0);
		tex_vec_vol->Bind(1);
		shader->Set_Uniform("tex3d",1);
		GLint tex_idx[4]={0,0,0,0};
		for(int i=0;i<d;i++){int idx=2+i;tex_trsf[i]->Bind(idx);tex_idx[i]=idx;}
		shader->Set_Uniform_Array("tex1d",4,&tex_idx[0]);
		glBindVertexArray(vao);
		glEnable(GL_CULL_FACE);
		glCullFace(GL_BACK);
		glEnable(GL_PRIMITIVE_RESTART);
		glPrimitiveRestartIndex((GLuint)UINT_MAX);
		glDrawElements(GL_TRIANGLE_STRIP,ele_size,GL_UNSIGNED_INT,0);
		glDisable(GL_PRIMITIVE_RESTART);
		glDisable(GL_CULL_FACE);
		Disable_Alpha_Blend();
		shader->End();}
    }

	virtual void Refresh(const int frame)
	{
		if(!initialized)Initialize();

		if(frame==0){std::string file_name=Object_File_Name(output_dir,frame,name);
			if(File::File_Exists(file_name)){File::Read_Binary_From_File(file_name,grid);}}
		
		Field<Vector<T,d>,3> field;
		for(int i=0;i<d;i++){std::string file_name=Object_File_Name(output_dir,frame,data[i].name);
			if(File::File_Exists(file_name)){
				std::shared_ptr<Field<real,3> > read_field=std::make_shared<Field<real,3> >();
				File::Read_Binary_From_File(file_name,*read_field);
				if(i==0){field.Resize(read_field->counts);field.Fill(Vector<T,d>::Zero());}

				iterate_cell_d(iter,grid,3){const Vector3i& cell=iter.Coord();
					int idx=field.counts[0]*field.counts[1]*cell[2]+field.counts[0]*cell[1]+cell[0];
					T val=(T)((float)0xffff*(float)(*read_field)(cell));field.array[idx][i]=val;}
				Set_Data_Refreshed();}}

		ArrayF<Array<T>,d> transfer;for(int i=0;i<d;i++)transfer[i]=c_trsf[i].transfer;

		Set_Data_Pointers(&grid,&field,&transfer);
	}
};
#endif