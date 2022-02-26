//////////////////////////////////////////////////////////////////////////
// Opengl uniform buffer object
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <GL/glew.h>
#include <GL/freeglut.h>
#include "MeshFunc.h"
#include "OpenGLUbos.h"
#include "OpenGLTexture.h"
#include "OpenGLPrimitives.h"
#include "OpenGLShaderLibrary.h"
#include "OpenGLMarkerObjects.h"

using namespace OpenGLUbos;
using namespace OpenGLTextures;

//////////////////////////////////////////////////////////////////////////
////OpenGLAxes

void OpenGLAxes::Initialize()	////No data update
{
	Base::Initialize();
	Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vcolor"));
	Update_Data_To_Render_Pre();
	ArrayF<OpenGLColor,3> colors={OpenGLColor::Red(),OpenGLColor::Green(),OpenGLColor::Blue()};
	int dim=use_2d_display?2:3;
	for(int i=0;i<dim;i++)for(int j=0;j<2;j++){
		Vector3 pos=Vector3::Zero()+Vector3::Unit(i)*axis_length*(real)j;
		OpenGL_Vertex4_And_Color4(pos,colors[i].rgba,opengl_vertices);}		////position, 4 floats; color, 4 floats
	Set_OpenGL_Vertices();
	Set_OpenGL_Vertex_Attribute(0,4,8,0);	////position
	Set_OpenGL_Vertex_Attribute(1,4,8,4);	////color
	Update_Data_To_Render_Post();
}

void OpenGLAxes::Display() const
{
	using namespace OpenGLUbos;
	if(!visible)return;
	Update_Polygon_Mode();
	{std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
	shader->Begin();
	Bind_Uniform_Block_To_Ubo(shader,"camera");
	glBindVertexArray(vao);
	glDrawArrays(GL_LINES,0,vtx_size/8);
	shader->End();}
}

//////////////////////////////////////////////////////////////////////////
////OpenGLPoint

void OpenGLPoint::Initialize()
{
	Base::Initialize();
	Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("psize_ucolor"));
}

void OpenGLPoint::Update_Data_To_Render()
{
	if(!Update_Data_To_Render_Pre())return;
	OpenGL_Vertex4(pos,opengl_vertices);	////position, 4 floats
	Set_OpenGL_Vertices();
	Set_OpenGL_Vertex_Attribute(0,4,4,0);	////position
	Update_Data_To_Render_Post();	
}

void OpenGLPoint::Display() const
{
	std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
	shader->Begin();
	Bind_Uniform_Block_To_Ubo(shader,"camera");
	shader->Set_Uniform_Vec4f("color",color.rgba);
	glEnable(GL_PROGRAM_POINT_SIZE);
	shader->Set_Uniform("point_size",point_size);
	glBindVertexArray(vao);
	glDrawArrays(GL_POINTS,0,vtx_size/4);
	shader->End();
}
	
//////////////////////////////////////////////////////////////////////////
////OpenGLTriangle

void OpenGLTriangle::Initialize()
{
	Base::Initialize();
	Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vpos"));
}

void OpenGLTriangle::Update_Data_To_Render()
{
	if(!Update_Data_To_Render_Pre())return;
	for(int i=0;i<3;i++)OpenGL_Vertex4(vtx[i],opengl_vertices);	////position, 4 floats
	Set_OpenGL_Vertices();
	Set_OpenGL_Vertex_Attribute(0,4,4,0);	////position
	Update_Data_To_Render_Post();	
}

void OpenGLTriangle::Display() const
{
	std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
	shader->Begin();
	glPushAttrib(GL_LINE_BIT);
	Bind_Uniform_Block_To_Ubo(shader,"camera");
	glLineWidth(line_width);
	shader->Set_Uniform_Vec4f("color",color.rgba);
	glBindVertexArray(vao);
	glDrawArrays(GL_LINE_LOOP,0,vtx_size/4);
	glPopAttrib();
	shader->End();
}

//////////////////////////////////////////////////////////////////////////
////OpenGLMarkerTriangleMesh

void OpenGLMarkerTriangleMesh::Initialize()
{
	Base::Initialize();
	Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vpos_model_vnormal_dl_fast"));

	Update_Mesh_Data_To_Render();
}

void OpenGLMarkerTriangleMesh::Update_Mesh_Data_To_Render()
{
	if(mesh.Vertices().size()==0)return;

	for(auto i=0;i<mesh.Vertices().size();i++){
		OpenGL_Vertex4(mesh.Vertices()[i],opengl_vertices);	////position, 4 floats
		OpenGL_Vertex4((*mesh.normals)[i],opengl_vertices);}
	Set_OpenGL_Vertices();
	Set_OpenGL_Vertex_Attribute(0,4,8,0);	////position
	Set_OpenGL_Vertex_Attribute(1,4,8,4);	////normal
		
	for(auto& e:mesh.elements)
		OpenGL_Vertex(e,opengl_elements);
	Set_OpenGL_Elements();
}

void OpenGLMarkerTriangleMesh::Update_Data_To_Render()
{
	if(!Update_Data_To_Render_Pre())return;
	Update_Model_Matrix();
	Update_Data_To_Render_Post();	
}

void OpenGLMarkerTriangleMesh::Display() const
{
	std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
	shader->Begin();
	Bind_Uniform_Block_To_Ubo(shader,"camera");
	Bind_Uniform_Block_To_Ubo(shader,"lights");
	shader->Set_Uniform_Vec4f("mat_dif",color.rgba);
	shader->Set_Uniform_Vec4f("mat_spec",color.rgba);
	shader->Set_Uniform_Matrix4f("model",glm::value_ptr(model));
	glBindVertexArray(vao);
	glDrawElements(GL_TRIANGLES,ele_size,GL_UNSIGNED_INT,0);
	shader->End();
}

//////////////////////////////////////////////////////////////////////////
////OpenGLSphere

void OpenGLSphere::Initialize()
{
	Base::Initialize();
	Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vpos_model_vnormal_dl_fast"));
	MeshFunc::Initialize_Sphere_Mesh(radius,&mesh);
	MeshFunc::Update_Normals(mesh);
	Update_Model_Matrix();
	Update_Mesh_Data_To_Render();
}

void OpenGLSphere::Update_Model_Matrix_Helper(const Vector3& pos,glm::mat4& model_matrix) const
{
	model_matrix=glm::translate(glm::mat4(1.f),glm::vec3((GLfloat)pos[0],(GLfloat)pos[1],(GLfloat)pos[2]));
}

void OpenGLSphere::Update_Model_Matrix()
{
	Update_Model_Matrix_Helper(pos,model);
}
	
void OpenGLSphere::Display_Multiple_Instances(const Array<Vector3>& centers) const
{
	std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
	shader->Begin();
	Bind_Uniform_Block_To_Ubo(shader,"camera");
	Bind_Uniform_Block_To_Ubo(shader,"lights");
	shader->Set_Uniform_Vec4f("mat_dif",color.rgba);
	shader->Set_Uniform_Vec4f("mat_spec",color.rgba);

	for(const auto& p:centers){
		glm::mat4 model_p;Update_Model_Matrix_Helper(p,model_p);
		shader->Set_Uniform_Matrix4f("model",glm::value_ptr(model_p));
		glBindVertexArray(vao);
		glDrawElements(GL_TRIANGLES,ele_size,GL_UNSIGNED_INT,0);}

	shader->End();
}

//////////////////////////////////////////////////////////////////////////
////OpenGLArrow

void OpenGLArrow::Initialize()
{
	Base::Initialize();
	Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vpos_model_vnormal_dl_fast"));
		
	real length=(real)1;
	real cylinder_r=std::min(length*(real).01,(real).01);
	real cylinder_h=length*(real).8;
	real cone_r=cylinder_r*(real)2;
	real cone_h=length-cylinder_h;
	int n=8;
	TriangleMesh<3> cylinder;MeshFunc::Initialize_Cylinder_Mesh(cylinder_r,cylinder_h,n,&cylinder);
	TriangleMesh<3> cone;MeshFunc::Initialize_Cone_Mesh(cone_r,cone_h,n,&cone);
	MeshFunc::Translate<3>(*cone.vertices,Vector3::Unit(2)*cylinder_h);

	Array<TriangleMesh<3>*> components;
	components.push_back(&cylinder);components.push_back(&cone);
	MeshFunc::Merge(components,&mesh);
	MeshFunc::Update_Normals(mesh);

	Update_Mesh_Data_To_Render();
}

void OpenGLArrow::Update_Model_Matrix()
{
	real angle=(real)0;Vector3 axis=Vector3::Unit(2);
	real s=(end-start).norm();
	AuxFunc::Angle_And_Axis_Between(Vector3::Unit(2),end-start,angle,axis);
	glm::mat4 rot=glm::rotate(glm::mat4(1.f),(GLfloat)angle,glm::vec3((GLfloat)axis[0],(GLfloat)axis[1],(GLfloat)axis[2]));
	glm::mat4 trans=glm::translate(glm::mat4(1.f),glm::vec3((GLfloat)start[0],(GLfloat)start[1],(GLfloat)start[2]));
	glm::mat4 scale=glm::scale(glm::mat4(1.f),glm::vec3(s,s,s));
	model=trans*rot*scale;
}

//////////////////////////////////////////////////////////////////////////
////OpenGLCircle

void OpenGLCircle::Initialize()
{
	Base::Initialize();
	Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vpos_model"));

	real step=two_pi/(real)n;
	for(int i=0;i<n;i++)OpenGL_Vertex4(Vector3(cos((real)i*step),sin((real)i*step),(real)0),opengl_vertices);	////position, 4 floats
	Set_OpenGL_Vertices();
	Set_OpenGL_Vertex_Attribute(0,4,4,0);	////position
}

void OpenGLCircle::Update_Data_To_Render()
{
	if(!Update_Data_To_Render_Pre())return;
	Update_Model_Matrix();
	Update_Data_To_Render_Post();
}

void OpenGLCircle::Display() const
{
	std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
	shader->Begin();
	glPushAttrib(GL_LINE_BIT);
	Bind_Uniform_Block_To_Ubo(shader,"camera");
	glLineWidth(line_width);
	shader->Set_Uniform_Vec4f("color",color.rgba);
	shader->Set_Uniform_Matrix4f("model",glm::value_ptr(model));
	glBindVertexArray(vao);
	switch(polygon_mode){
	case PolygonMode::Fill:
		glDrawArrays(GL_POLYGON,0,vtx_size/4);break;
	case PolygonMode::Wireframe:
		glDrawArrays(GL_LINE_LOOP,0,vtx_size/4);break;}
	glPopAttrib();
	shader->End();
}

void OpenGLCircle::Display_Multiple_Instances(const Array<Vector3>& centers) const
{
	std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
	shader->Begin();
	glPushAttrib(GL_LINE_BIT);
	Bind_Uniform_Block_To_Ubo(shader,"camera");
	glLineWidth(line_width);
	shader->Set_Uniform_Vec4f("color",color.rgba);

	for(const auto& p:centers){
		glm::mat4 model_p;Update_Model_Matrix_Helper(p,radius,model_p);
		shader->Set_Uniform_Matrix4f("model",glm::value_ptr(model_p));
		glBindVertexArray(vao);
		switch(polygon_mode){
		case PolygonMode::Fill:
			glDrawArrays(GL_POLYGON,0,vtx_size/4);break;
		case PolygonMode::Wireframe:
			glDrawArrays(GL_LINE_LOOP,0,vtx_size/4);break;}}

	glPopAttrib();
	shader->End();
}

void OpenGLCircle::Display_Multiple_Instances(const Array<Vector3>& centers,const Array<real>& radii) const
{
	std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
	shader->Begin();
	glPushAttrib(GL_LINE_BIT);
	Bind_Uniform_Block_To_Ubo(shader,"camera");
	glLineWidth(line_width);
	shader->Set_Uniform_Vec4f("color",color.rgba);

	for(int i=0;i<(int)centers.size();i++){
		glm::mat4 model_p;Update_Model_Matrix_Helper(centers[i],radii[i],model_p);
		shader->Set_Uniform_Matrix4f("model",glm::value_ptr(model_p));
		glBindVertexArray(vao);
		switch(polygon_mode){
		case PolygonMode::Fill:
			glDrawArrays(GL_POLYGON,0,vtx_size/4);break;
		case PolygonMode::Wireframe:
			glDrawArrays(GL_LINE_LOOP,0,vtx_size/4);break;}}

	glPopAttrib();
	shader->End();
}

void OpenGLCircle::Update_Model_Matrix()
{
	Update_Model_Matrix_Helper(pos,radius,model);
}

void OpenGLCircle::Update_Model_Matrix_Helper(const Vector3& pos,const real radius,glm::mat4& model_matrix) const
{
	model_matrix=glm::translate(glm::mat4(1.f),glm::vec3((GLfloat)pos[0],(GLfloat)pos[1],(GLfloat)pos[2]));
	model_matrix=glm::scale(model_matrix,glm::vec3(radius,radius,radius));
}

void OpenGLSquare::Initialize()
{
	Base::Initialize();
	Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vpos_model"));

	OpenGL_Vertex4(Vector3(0,0,(real)0),opengl_vertices);	////position, 4 floats
	OpenGL_Vertex4(Vector3(axis_length,0,(real)0),opengl_vertices);	////position, 4 floats
	OpenGL_Vertex4(Vector3(axis_length,axis_length,(real)0),opengl_vertices);	////position, 4 floats
	OpenGL_Vertex4(Vector3(0,axis_length,(real)0),opengl_vertices);	////position, 4 floats

	Set_OpenGL_Vertices();
	Set_OpenGL_Vertex_Attribute(0,4,4,0);	////position
	Update_Model_Matrix();
}

void OpenGLSquare::Display() const
{
	std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
	shader->Begin();
	glPushAttrib(GL_LINE_BIT);
	Bind_Uniform_Block_To_Ubo(shader,"camera");
	glLineWidth(line_width);
	shader->Set_Uniform_Vec4f("color",color.rgba);
	shader->Set_Uniform_Matrix4f("model",glm::value_ptr(model));
	glBindVertexArray(vao);
	glDrawArrays(GL_LINE_LOOP,0,vtx_size/4);
	glPopAttrib();
	shader->End();
}

void OpenGLSquare::Update_Model_Matrix()
{
	model=glm::translate(glm::mat4(1.f),glm::vec3((GLfloat)pos[0],(GLfloat)pos[1],(GLfloat)pos[2]));
}

void OpenGLSegments::Initialize()
{
	Base::Initialize();
	Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vpos_model"));

	for(int i=0;i<vertices.size();i++)
		OpenGL_Vertex4(vertices[i],opengl_vertices);	////position, 4 floats

	Set_OpenGL_Vertices();
	Set_OpenGL_Vertex_Attribute(0,4,4,0);	////position
}

void OpenGLSegments::Display() const
{
	std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
	shader->Begin();
	glPushAttrib(GL_LINE_BIT);
	Bind_Uniform_Block_To_Ubo(shader,"camera");
	glLineWidth(line_width);
	shader->Set_Uniform_Vec4f("color",color.rgba);
	shader->Set_Uniform_Matrix4f("model",glm::value_ptr(model));
	glBindVertexArray(vao);
	glDrawArrays(GL_LINES,0,vtx_size/4);
	glPopAttrib();
	shader->End();
}