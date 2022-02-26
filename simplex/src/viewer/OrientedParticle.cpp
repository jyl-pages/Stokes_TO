#include "OrientedParticle.h"

void Render_Tracker_Particles(const Array<Vector3>& pos, std::shared_ptr<OpenGLShaderProgram> shader, const OpenGLColor &color, GLfloat point_size) {
	static Array<GLfloat> vertices;
	vertices.resize(pos.size() * 4);
	for (int i = 0; i < pos.size(); i++) {
		const Vector3& p = pos[i];
		vertices[i * 4] = (GLfloat)p[0];
		vertices[i * 4 + 1] = (GLfloat)p[1];
		vertices[i * 4 + 2] = (GLfloat)p[2];
		vertices[i * 4 + 3] = (GLfloat)0;
	}
	GLuint vao = 0;
	using namespace OpenGLUbos; using namespace OpenGLTextures;
	shader->Begin();
	Bind_Uniform_Block_To_Ubo(shader, "camera");
	shader->Set_Uniform_Vec4f("color", color.rgba);
	glEnable(GL_PROGRAM_POINT_SIZE);
	shader->Set_Uniform("point_size", point_size);
	glBindVertexArray(vao);
	glDrawArrays(GL_POINTS, 0, pos.size());
	shader->End();
}

void OrientedParticle::Initialize(void)
{
	Base::Initialize();
	Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("psize_ucolor"));
	Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vpos_model"));	//In OpenGLShaderLibrary.cpp
	Add_Shader_Program(OpenGLShaderLibrary::Get_Shader_From_File("TrackerCircles", Path::Data() + "/shaders/circle.vert", Path::Data() + "/shaders/circle.frag"));
}

void OrientedParticle::Initialize_Data(std::string _output_dir, const int frame)
{
	output_dir = _output_dir;
	Refresh(frame);
}

void OrientedParticle::Refresh(const int frame)
{
	std::cout << "refresh " << frame << std::endl;
	std::string frame_path = output_dir + "/" + std::to_string(frame);
	pos.Read_Data(frame_path);
	normal.Read_Data(frame_path);
	for (int i = 0; i < scalar_qnts.size(); i++) scalar_qnts[i].Read_Data(frame_path);
	for (int i = 0; i < vector_qnts.size(); i++) vector_qnts[i].Read_Data(frame_path);
	Set_Data_Refreshed();
}

void OrientedParticle::Update_Data_To_Render(void)
{
	if (!Update_Data_To_Render_Pre())return;
	mat.mat_dif = glm::vec4(color.rgba[0], color.rgba[1], color.rgba[2], color.rgba[3]);
	Set_Data_Refreshed(false);
}

void OrientedParticle::Update_Model_Matrix_Helper(const float* pos, const float* dir, const float radius, glm::mat4& model_matrix) const
{
	model_matrix = glm::translate(glm::mat4(1.f), glm::vec3(pos[0], pos[1], pos[2]));
	float rad; Vector3f axis; AuxFunc::Angle_And_Axis_Between(Vector3f(0.f, 0.f, 1.f), Vector3f(dir[0], dir[1], dir[2]), rad, axis);
	model_matrix = glm::rotate(model_matrix, rad, glm::vec3(axis[0], axis[1], axis[2]));
	model_matrix = glm::scale(model_matrix, glm::vec3(radius, radius, radius));
}

void OrientedParticle::Display_Tracker_Points(void)
{
	static Array<GLfloat> vertices;
	vertices.resize(pos.data.size() * 4);
	for (int i = 0; i < pos.data.size(); i++) {
		const Vector3& p = pos.data[i];
		vertices[i * 4] = (GLfloat)p[0];
		vertices[i * 4 + 1] = (GLfloat)p[1];
		vertices[i * 4 + 2] = (GLfloat)p[2];
		vertices[i * 4 + 3] = 0;
	}

}

void OrientedParticle::Display() const
{
	Update_Polygon_Mode();
	Render_Tracker_Particles(pos.data, shader_programs[0], pos.color, point_size);
}

void OrientedParticle::Add_Scalar(std::string _name, bool _visible, OpenGLColor _color)
{
	LagrangianScalar s;
	s.Initialize(_name, _visible, _color);
	scalar_qnts.push_back(s);
}

void OrientedParticle::Add_Vector(std::string _name, bool _visible, OpenGLColor _color)
{
	LagrangianVector v;
	v.Initialize(_name, _visible, _color);
	vector_qnts.push_back(v);
}


void LagrangianQuantity::Initialize(std::string _filename, bool _visible, OpenGLColor _color)
{
	file_name = _filename;
	visible = _visible;
	color = _color;
}

void LagrangianScalar::Read_Data(const std::string& frame_path)
{
	BinaryDataIO::Read_Scalar_Array(frame_path + "/" + file_name, data);
}

void LagrangianVector::Read_Data(const std::string& frame_path)
{
	BinaryDataIO::Read_Vector_Array_3D<real, 3>(frame_path + "/" + file_name, data);
}

