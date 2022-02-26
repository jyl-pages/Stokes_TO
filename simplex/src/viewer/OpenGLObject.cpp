//////////////////////////////////////////////////////////////////////////
// Opengl object
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "AuxFunc.h"
#include "File.h"
#include "OpenGLShaderProgram.h"
#include "OpenGLObject.h"
#include "OpenGLPrimitives.h"
#include "OpenGLColorMapper.h"

OpenGLObject::OpenGLObject()
{
	color_mapper.reset(new OpenGLColorMapper());
	Initialize_Callbacks();
}

void OpenGLObject::Initialize_OpenGL_Buffers()
{
	glGenVertexArrays(1,&vao);
	glGenBuffers(1,&vbo);
	glGenBuffers(1,&ebo);
}

void OpenGLObject::Set_Env_Mapping(const std::string& _env_name)
{
	env_name=_env_name;
	Set_Shading_Mode(ShadingMode::EnvMapping);
}

bool OpenGLObject::Update_Data_To_Render_Pre()
{
	if(!initialized)Initialize();
	if(!data_refreshed)return false;
	Update_Color_Mapper();
	Clear_OpenGL_Arrays();
	return true;
}

void OpenGLObject::Update_Data_To_Render_Post()
{
	Clear_OpenGL_Arrays();
	Set_Data_Refreshed(false);
}

////io helper functions
std::string OpenGLObject::Object_File_Name(const std::string& output_dir,const int frame,const std::string& object_name)
{return output_dir+"/"+std::to_string(frame)+"/"+object_name;}

bool OpenGLObject::Object_File_Exists(const std::string& output_dir,const int frame,const std::string& object_name)
{std::string file_name=Object_File_Name(output_dir,frame,object_name);return File::File_Exists(file_name);}

void OpenGLObject::Output_Debug_Info(std::ostream & out)
{
	out << name << ", visible: " << visible;
}

////Opengl helper functions
void OpenGLObject::Update_Polygon_Mode() const
{
	switch(polygon_mode){
	case PolygonMode::Fill:glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);break;
	case PolygonMode::Wireframe:
	case PolygonMode::SurfOnly:{glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);}break;}
}

void OpenGLObject::Clear_OpenGL_Arrays()
{
	opengl_vertices.clear();
	opengl_elements.clear();
}

void OpenGLObject::Set_OpenGL_Vertices()
{
	vtx_size=(int)opengl_vertices.size(); 
	Set_OpenGL_Vertices(opengl_vertices,vtx_size);
}

void OpenGLObject::Set_OpenGL_Vertices(const Array<GLfloat>& opengl_vertices,int vtx_size,int offset)const
{
	if (opengl_vertices.empty()) return;//an empty Array does not have a [0]
	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER,vbo);
	glBufferData(GL_ARRAY_BUFFER,opengl_vertices.size()*sizeof(GLfloat),&opengl_vertices[offset],GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER,0);
	glBindVertexArray(0);	
}

void OpenGLObject::Set_OpenGL_Elements()
{Set_OpenGL_Elements(opengl_elements,ele_size);}

void OpenGLObject::Set_OpenGL_Elements(Array<GLuint>& opengl_elemetns,int& ele_size)
{
	glBindVertexArray(vao);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,ebo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,opengl_elements.size()*sizeof(GLuint),&opengl_elements[0],GL_STATIC_DRAW);
	glBindVertexArray(0);
	ele_size=(int)opengl_elements.size();
}

void OpenGLObject::Set_OpenGL_Vertex_Attribute(const GLuint idx,const GLint element_size,const GLuint stride_size,GLuint start_idx)
{
	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER,vbo);
	glVertexAttribPointer(idx,element_size,GL_FLOAT,GL_FALSE,stride_size*sizeof(GLfloat),(GLvoid*)(start_idx*sizeof(GLfloat)));
	glEnableVertexAttribArray(idx);
	glBindBuffer(GL_ARRAY_BUFFER,0);
	glBindVertexArray(0);	
}

void OpenGLObject::Set_Material_Color(const OpenGLColor& c)
{
	use_mat=true;
	mat.mat_dif=glm::vec4(c.rgba[0],c.rgba[1],c.rgba[2],c.rgba[3]);
	mat.mat_spec=glm::vec4(c.rgba[0],c.rgba[1],c.rgba[2],c.rgba[3]);
}

bool OpenGLObject::Use_Alpha_Blend() const {return alpha<1.f;}

void OpenGLObject::Enable_Alpha_Blend() const
{glEnable(GL_BLEND);glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);glBlendEquation(GL_FUNC_ADD);}

void OpenGLObject::Disable_Alpha_Blend() const {glDisable(GL_BLEND);}

void OpenGLObject::Update_Scalar_Range(const Array<Matrix3>& array,real& v_min,real& v_max) const
{
	v_min=std::numeric_limits<real>::max();v_max=std::numeric_limits<real>::min();
	for(auto& v:array){Sym_Mat3_Eig eig(v,false);
		for(int i=0;i<3;i++){real s=abs(eig.eigenvalues()[i]);
			if(s<v_min)v_min=s;else if(s>v_max)v_max=s;}}
}

real OpenGLObject::Scalar(const real& v) const {return v;}
real OpenGLObject::Scalar(const Vector3& v) const {return v.norm();}
real OpenGLObject::Scalar(const Matrix3& v) const {Sym_Mat3_Eig eig(v,false);return AuxFunc::Abs_Max(eig.eigenvalues());}

void OpenGLObject::Toggle_Draw()
{visible=!visible;if(verbose)std::cout<<"Toggle visible: "<<visible<<", ["<<name<<"]"<<std::endl;}

void OpenGLObject::Toggle_Normalize()
{normalize=!normalize;if(verbose)std::cout<<"Toggle normalize: "<<normalize<<", ["<<name<<"]"<<std::endl;}

void OpenGLObject::Toggle_Draw_Dis()
{draw_dis=!draw_dis;if(verbose)std::cout<<"Toggle draw_dis: "<<draw_dis<<", ["<<name<<"]"<<std::endl;}

void OpenGLObject::Toggle_Increase_Scale()
{real rescale=(real)1.1;scale*=rescale;Set_Data_Refreshed();}

void OpenGLObject::Toggle_Decrease_Scale()
{real rescale=(real)1.1;scale/=rescale;Set_Data_Refreshed();}


void OpenGLObject::Toggle_Increase_Feature(){}

void OpenGLObject::Toggle_Decrease_Feature(){}

void OpenGLObject::Initialize_Callbacks()
{
	callbacks.push_back(&Toggle_Data_Idx_0_Func);
	callbacks.push_back(&Toggle_Data_Idx_1_Func);
	callbacks.push_back(&Toggle_Data_Idx_2_Func);
	callbacks.push_back(&Toggle_Data_Idx_3_Func);

	callbacks.push_back(&Toggle_Func_Idx_0_Func);
	callbacks.push_back(&Toggle_Func_Idx_1_Func);
	callbacks.push_back(&Toggle_Func_Idx_2_Func);
	callbacks.push_back(&Toggle_Func_Idx_3_Func);
}