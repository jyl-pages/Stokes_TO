//////////////////////////////////////////////////////////////////////////
// Opengl uniform buffer object
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <GL/glew.h>
#include <GL/freeglut.h>
#include "GeometryPrimitives.h"
#include "Mesh.h"
#include "MeshFunc.h"
#include "OpenGLUbos.h"
#include "OpenGLTexture.h"
#include "OpenGLPrimitives.h"
#include "OpenGLShaderLibrary.h"
#include "OpenGLColorMapper.h"
#include "OpenGLScreenObjects.h"
#include "OpenGLWindow.h"

//////////////////////////////////////////////////////////////////////////
////OpenGLScreenObject
using namespace OpenGLUbos;
using namespace OpenGLTextures;

void OpenGLScreenObject::Initialize()
{
	Base::Initialize();
	Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("ortho_vcolor"));
}

void OpenGLScreenObject::Set_Pos(const Vector2& _pos)
{pos=_pos;model=glm::translate(glm::mat4(),glm::vec3((GLfloat)pos[0],(GLfloat)pos[1],0.f));}

void OpenGLText::Display() const
{
	if(visible==false||texts.size()==0||ortho==nullptr)return;

	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadMatrixf(glm::value_ptr(*ortho));
	glColor3f(color.rgba[0],color.rgba[1],color.rgba[2]);
		
	for(auto i=0;i<texts.size();i++){
		Vector2 start=pos+offsets[i];
		glRasterPos2d(start[0],start[1]);
		const std::string& str=texts[i];
		for(auto j=0;j<str.length();j++) glutBitmapCharacter(GLUT_BITMAP_9_BY_15,str[j]);}

	glPopAttrib();
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

Vector2 OpenGLText3D::Proj_To_Win_Pos(const Vector3& pos_3d) const
{
	Vector2 pos=Vector2::Zero();
	if(opengl_window!=nullptr){
		Vector3f win_bl_pos=opengl_window->Project(pos_3d.cast<float>());
		pos=Vector2((real)win_bl_pos[0],(real)win_bl_pos[1]);}
	return pos;
}

void OpenGLText3D::Set_Pos(const Vector3& _pos_3d)
{
	pos_3d=_pos_3d;
}

void OpenGLText3D::Display() const
{
	if(visible==false||texts.size()==0||ortho==nullptr)return;
	Vector2 pos_2d=Proj_To_Win_Pos(pos_3d);

	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadMatrixf(glm::value_ptr(*ortho));
	glColor3f(color.rgba[0],color.rgba[1],color.rgba[2]);
		
	for(auto i=0;i<texts.size();i++){
		Vector2 start=pos_2d+offsets[i];
		glRasterPos2d(start[0],start[1]);
		const std::string& str=texts[i];
		for(auto j=0;j<str.length();j++) glutBitmapCharacter(GLUT_BITMAP_9_BY_15,str[j]);}

	glPopAttrib();
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();	
}

void OpenGLTextArray3D::Display() const
{
	if(visible==false||texts.size()==0||ortho==nullptr)return;

	glPushAttrib(GL_ENABLE_BIT);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadMatrixf(glm::value_ptr(*ortho));
	glColor3f(color.rgba[0],color.rgba[1],color.rgba[2]);
		
	for(auto i=0;i<texts.size();i++){
		Vector2 pos_2d=Proj_To_Win_Pos(pos_3d_array[i]);
		Vector2 start=pos_2d+offsets[i];
		glRasterPos2d(start[0],start[1]);
		const std::string& str=texts[i];
		for(auto j=0;j<str.length();j++) glutBitmapCharacter(GLUT_BITMAP_9_BY_15,str[j]);}

	glPopAttrib();
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();	
}

void OpenGLTextArray3D::Update_Vertex_Data(const Array<Vector3>& vertices)
{
	int n=(int)vertices.size();pos_3d_array.resize(n);texts.resize(n);offsets.resize(n);
	for(int i=0;i<n;i++){texts[i]=std::to_string(i);pos_3d_array[i]=vertices[i];offsets[i]=Vector2::Zero();}
}

void OpenGLTextArray3D::Update_Vertex_Data(const Array<Vector3>& vertices,const Hashset<int>& selected_vertices)
{
	int n=(int)vertices.size();pos_3d_array.resize(n);texts.resize(n);offsets.resize(n);
	for(int i=0;i<n;i++){if(selected_vertices.find(i)==selected_vertices.end())continue;
		texts.push_back(std::to_string(i));pos_3d_array.push_back(vertices[i]);offsets.push_back(Vector2::Zero());}
}

void OpenGLTextArray3D::Update_Element_Data(const Array<Vector3>& vertices,const Array<Vector3i>& elements)
{
	int n=(int)elements.size();pos_3d_array.resize(n);texts.resize(n);offsets.resize(n);
	for(int i=0;i<n;i++){texts[i]=std::to_string(i);pos_3d_array[i]=MeshFunc::Element_Center(vertices,elements,i);offsets[i]=Vector2::Zero();}
}

//////////////////////////////////////////////////////////////////////////
////OpenGLBar

void OpenGLBar::Initialize()
{
	Base::Initialize();
	Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("ortho_vcolor"));
	color_mapper->Reset_Color_Mode((real)0,(real)1,0);
}

void OpenGLBar::Update_Data_To_Render()
{
	if(!Update_Data_To_Render_Pre())return;
	Set_Bar_Mesh();

	Set_OpenGL_Vertices();
	Set_OpenGL_Vertex_Attribute(0,2,8,0);	////position
	Set_OpenGL_Vertex_Attribute(1,4,8,4);	////color
	Set_OpenGL_Elements();
	Update_Data_To_Render_Post();
}

void OpenGLBar::Set_Bar_Mesh()
{
	real step=bar_width/(real)res;
	Vector2 left_top=box.min_corner+Vector2::Unit(1)*bar_height;
	int e=0;for(int i=0;i<=res;i++)for(int j=0;j<2;j++){
		Vector2 pos=left_top+Vector2::Unit(0)*(real)i*step-Vector2::Unit(1)*(real)j*bar_height;
		//real ratio=(pos[0]-bar.min_corner[0])/width;
		//OpenGLColor c=color_mapper->Color(ratio);
		OpenGL_Vertex4_And_Color4(pos,color.rgba,opengl_vertices);		////position, 4 floats; color, 4 floats
		OpenGL_Vertex(e++,opengl_elements);}	
}

void OpenGLBar::Display() const
{
	if(!visible)return;

	Update_Polygon_Mode();
	{std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
	shader->Begin();
	OpenGLUbos::Bind_Uniform_Block_To_Ubo(shader,"camera");
	glBindVertexArray(vao);
	shader->Set_Uniform_Matrix4f("model",glm::value_ptr(model));
	glDrawArrays(GL_TRIANGLE_STRIP,0,vtx_size/8);
	shader->End();}
}

//////////////////////////////////////////////////////////////////////////
////OpenGLColorBar

void OpenGLColorBar::Initialize()
{
	Base::Initialize();

	glm::mat4* ortho=nullptr;Camera* camera=Get_Camera();
	if(camera!=nullptr){ortho=&camera->ortho;}
	scalar_texts.Set_Data_Pointers(ortho);
}

void OpenGLColorBar::Set_Data_Pointer(OpenGLColorMapper* _color_mapper_ptr,int* _data_idx_ptr)
{color_mapper_ptr=_color_mapper_ptr;data_idx_ptr=_data_idx_ptr;}

void OpenGLColorBar::Set_Bar_Mesh()
{
	real step=bar_height/(real)res;
	Vector2 left_bot=box.min_corner;
	real scale_width=bar_width*(real).5;

	if(draw_scalars){
		scalar_vertices.clear();
		scalar_texts.Clear_Data();}

	int e=0;for(int i=0;i<=res;i++)for(int j=0;j<2;j++){
		real x=(real)j*bar_width;real y=(real)i*step;
		Vector2 pos=left_bot+Vector2(x,y);
			
		real ratio=(pos[1]-left_bot[0])/bar_height;
		OpenGLColor c;real scalar_value=ratio;
		if(color_mapper_ptr!=nullptr&&data_idx_ptr!=nullptr){
			scalar_value=ratio*(color_mapper_ptr->range[color_mapper_ptr->current_mode_index][1]-
				color_mapper_ptr->range[color_mapper_ptr->current_mode_index][0])
				+color_mapper_ptr->range[color_mapper_ptr->current_mode_index][0];
			c=color_mapper_ptr->Color(scalar_value);}
		else c=color_mapper->Color(ratio);
		//else c=color;

		OpenGL_Vertex4_And_Color4(pos,c.rgba,opengl_vertices);		////position, 4 floats; color, 4 floats
		OpenGL_Vertex(e++,opengl_elements);

		if(draw_scalars){
			real sx=bar_width+(real)j*scale_width;
			if((i==0||i==res)&&j==0)sx=(real)0;
			Vector2 scale_pos=left_bot+Vector2(sx,y);
			scalar_vertices.push_back(scale_pos);

			if(j==1){
				real space=bar_width;
				std::string str=std::to_string(scalar_value);
				str.erase(str.find_last_not_of('0')+1,std::string::npos);
				scalar_texts.texts.push_back(str);
				scalar_texts.offsets.push_back(Vector2(x+space,y));}}}

	bar_vtx_size=(int)opengl_vertices.size();
	for(int i=0;i<(int)scalar_vertices.size();i++){
		OpenGL_Vertex4(scalar_vertices[i],opengl_vertices);
		OpenGL_Color4(OpenGLColor::Gray(.5f).rgba,opengl_vertices);
		OpenGL_Vertex(e++,opengl_elements);}
	scale_vtx_size=(int)opengl_vertices.size()-bar_vtx_size;
}

void OpenGLColorBar::Display() const
{
	if(!visible)return;

	Update_Polygon_Mode();
	{std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
	shader->Begin();
	Bind_Uniform_Block_To_Ubo(shader,"camera");
	glBindVertexArray(vao);
	shader->Set_Uniform_Matrix4f("model",glm::value_ptr(model));
	glDrawArrays(GL_TRIANGLE_STRIP,0,bar_vtx_size/8);
	if(draw_scalars)glDrawArrays(GL_LINES,bar_vtx_size/8,scale_vtx_size/8);
	shader->End();}

	if(draw_scalars)scalar_texts.Display();
}

//////////////////////////////////////////////////////////////////////////
////OpenGLBackground

OpenGLBackground::OpenGLBackground()
{
	color=OpenGLColor::White();name="background";
	box=Box<2>(Vector2::Ones()*(real)-1,Vector2::Ones());polygon_mode=PolygonMode::Fill;
	//Set_Texture("ocean_blue.jpg");	////comment back to set up the texture background
	Set_Depth((real).9999);
}

void OpenGLBackground::Initialize()
{
	Base::Initialize();
	//mix_colors={OpenGLColor(.0f,.0f,.0f,1.f),OpenGLColor(.0f,.1f,.0f,1.f)};		////default background color: dark
	//mix_colors={OpenGLColor(1.f,1.f,1.f,1.f),OpenGLColor(.9f,1.f,.9f,1.f)};	////default background color: light
	mix_colors={OpenGLColor::Gray(.4f),OpenGLColor::Gray(.2f)};
	Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("gcolor_bk"));	////gradient color
	Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vtex_bk"));

	Update_Data_To_Render_Pre();

	Array<Vector3> vtx={Vector3(box.min_corner[0],box.min_corner[1],depth),
						Vector3(box.max_corner[0],box.min_corner[1],depth),
						Vector3(box.max_corner[0],box.max_corner[1],depth),
						Vector3(box.min_corner[0],box.max_corner[1],depth)};
	Array<Vector2> uv={Vector2((real)0.,(real)0.),Vector2((real)1.,(real)0.),Vector2((real)1.,(real)1.),Vector2((real)0.,(real)1.)};

	if(!use_vtx_tex)for(auto& p:vtx){
		OpenGL_Vertex4(p,opengl_vertices);			////position, 4 floats
		OpenGL_Color4(color.rgba,opengl_vertices);}	////color, 4 floats
	else for(auto i=0;i<vtx.size();i++){
		OpenGL_Vertex4(vtx[i],opengl_vertices);		////position, 4 floats
		OpenGL_Vertex4(uv[i],opengl_vertices);}		//// texture, 4 floats

	Set_OpenGL_Vertices();
	Set_OpenGL_Vertex_Attribute(0,4,8,0);	////position
	Set_OpenGL_Vertex_Attribute(1,4,8,4);	////color/tex
		
	Update_Data_To_Render_Post();	
}

void OpenGLBackground::Display() const
{
    using namespace OpenGLTextures;
	if(!visible)return;
	Update_Polygon_Mode();

	{std::shared_ptr<OpenGLShaderProgram> shader=use_vtx_tex?shader_programs[1]:shader_programs[0];
	shader->Begin();
	glDepthMask(GL_FALSE);
	if(use_vtx_tex){
		OpenGLTextures::Bind_Texture(tex_name,TextureType::Tx2d);
		shader->Set_Uniform("tex2d",0);}
	else{
		for(int i=0;i<mix_colors.size();i++){
			std::string name="color_"+std::to_string(i);
			shader->Set_Uniform_Vec4f(name,mix_colors[i].rgba);}}
	glBindVertexArray(vao);
	glDrawArrays(GL_QUADS,0,vtx_size/8);
	glDepthMask(GL_TRUE);
	shader->End();}
}