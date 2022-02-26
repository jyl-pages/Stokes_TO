//////////////////////////////////////////////////////////////////////////
// Opengl tracker vectors
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLTrackerCircles_h__
#define __OpenGLTrackerCircles_h__
#include <numeric>
#include "Common.h"
#include "Particles.h"
#include "File.h"
#include "AuxFunc.h"
#include "OpenGLObject.h"
#include "OpenGLShaderLibrary.h"
#include "OpenGLUbos.h"
#include "OpenGLTexture.h"
#include "OpenGLPrimitives.h"

#include "ArrayIO.h"

//////////////////////////////////////////////////////////////////////////
////This class is for the visualization of large-scale vector clouds

class OpenGLTrackerCircles : public OpenGLObject
{
	typedef OpenGLObject Base;
public:
	bool use_var_radius = false;
	float radius = (float).02;
	glm::mat4 model = glm::mat4(1.f);
	int n = 8;
	Array<GLfloat> circle_vertices;
	Array<GLfloat> point_vertices;

	Array<Vector3d> color_table;
	real max_color_thickness = -1.;
	real best_middle_color = 300.;
	bool first_time = true;
	real last_scale = 1.;
	real color_multiplier;
	std::string height_file_name = "h_bin";
	std::string sa_file_name = "";
	std::string color_table_file_name = "./color_table.dat";
	Array<real> sas;

	OpenGLTrackerCircles() { color = OpenGLColor::Red(); name = "tracker_circles"; }

	void Set_SA_File(std::string _name) { sa_file_name = _name; }
	void Set_Height_File(std::string _name) { height_file_name = _name; }
	void Set_Color_File(std::string _name) { color_table_file_name = _name; }

	void Load_Color_Table(void) {
		//load color table
		double color_r, color_g, color_b;
		std::ifstream file;
		file.open(color_table_file_name, std::ios::binary);
		if (!file.is_open()) {
			std::cerr << "[Warning]OpenGLTrackerColoredCircles: no valid color table file\n";
			return;
		}
		else {
			while (file.read(reinterpret_cast<char*>(&color_r), sizeof(double))
				&& file.read(reinterpret_cast<char*>(&color_g), sizeof(double))
				&& file.read(reinterpret_cast<char*>(&color_b), sizeof(double))) {
				Vector3d curr_color;
				curr_color << color_r, color_g, color_b;
				color_table.push_back(curr_color);
				max_color_thickness += 1.;
			}
			file.close();
			std::cout << "Color table goes up to: " << max_color_thickness << std::endl;
		}
	}

	virtual void Initialize()
	{
		Base::Initialize();
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vpos_model"));	//In OpenGLShaderLibrary.cpp
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader_From_File("TrackerCircles", Path::Data() + "/shaders/circle.vert", Path::Data() + "/shaders/circle.frag"));

		float step = (float)two_pi / (float)n;

		switch (shading_mode) {
		case ShadingMode::Lighting: {
			////default axis is z (0,0,1)
			for (int i = 0; i < n; i++) {
				OpenGL_Vertex4(Vector3(cos((float)i * step), sin((float)i * step), (float)0), circle_vertices);	////position, 4 floats
				OpenGL_Vertex4(Vector3(0., 0., 1.), circle_vertices);
			}											////normal, 4 floats
			vtx_size = circle_vertices.size();
			Base::Set_OpenGL_Vertices(circle_vertices, vtx_size);
			Set_OpenGL_Vertex_Attribute(0, 4, 8, 0);	////position
			Set_OpenGL_Vertex_Attribute(1, 4, 8, 4);	////normal	
		}break;
		case ShadingMode::Multicolor: {
			////default axis is z (0,0,1)
			for (int i = 0; i < n; i++) {
				OpenGL_Vertex4(Vector3(0.4 * cos((float)i * step), 0.4 * sin((float)i * step), 0.4 * (float)0), circle_vertices);	////position, 4 floats
				OpenGL_Vertex4(Vector3(0., 0., 1.), circle_vertices);
			}											////normal, 4 floats
			vtx_size = circle_vertices.size();
			Base::Set_OpenGL_Vertices(circle_vertices, vtx_size);
			Set_OpenGL_Vertex_Attribute(0, 4, 8, 0);	////position
			Set_OpenGL_Vertex_Attribute(1, 4, 8, 4);	////normal
		}break;
		default: {
			////default axis is z (0,0,1)
			for (int i = 0; i < n; i++) {
				OpenGL_Vertex4(Vector3(cos((float)i * step), sin((float)i * step), (float)0), circle_vertices);
			}	////position, 4 floats	
			Set_OpenGL_Vertex_Attribute(0, 4, 4, 0);	////position
		}break;
		}
	}

	virtual void Update_Data_To_Render()
	{
		if (!Update_Data_To_Render_Pre())return;
		mat.mat_dif = glm::vec4(color.rgba[0], color.rgba[1], color.rgba[2], color.rgba[3]);
		Update_Data_To_Render_Post();
	}

	virtual void Refresh(const int frame)
	{
		std::string file_name = output_dir + "/" + std::to_string(frame) + "/" + name;
		if (!File::File_Exists(file_name)) { return; }//just silently ignore me, do not crash
		std::ifstream input(file_name, std::ios::binary);
		if (!input) {
			std::cerr << "OpenGLTrackerCircles error: cannot open file " << file_name << "\n";
			exit(1);
		}
		int n = 0; File::Read_Binary(input, n);
		point_vertices.resize(n);
		File::Read_Binary_Array(input, &point_vertices[0], n);

		if (shading_mode == ShadingMode::Multicolor) {
			//read height
			Array<real> hs;
			std::string h_name = output_dir + "/" + std::to_string(frame) + "/" + height_file_name;
			if (!File::File_Exists(h_name)) { return; }
			BinaryDataIO::Read_Scalar_Array<real>(h_name, hs);
			int pn = hs.size();
			if (first_time) {
				if (max_color_thickness < 0) Load_Color_Table();
				if (max_color_thickness < 0) { std::cerr << "OpenGLTrackerCircles error: didn't load color table before\n"; exit(1); }
				real sum = 0;
				for (int i = 0; i < pn; i++) {
					sum += hs[i];
				}
				real avg_height = (real)(1. / hs.size()) * sum;
				if (avg_height < 1e-8) { avg_height = 1e-8; }
				if (color_multiplier == 0.) color_multiplier = best_middle_color / avg_height;
				if (sum != 0.)first_time = false;
			}
			//map height to colors
			multi_color.clear();
			for (int i = 0; i < pn; i++) {
				int curr_h = (int)(hs[i] * color_multiplier);
				Vector3d curr_rgb;
				if (curr_h <= max_color_thickness && curr_h >= 0) { curr_rgb = color_table[curr_h]; }
				else if (curr_h > max_color_thickness) { curr_rgb = Vector3d::Ones(); }
				else { curr_rgb = Vector3d::Zero(); }
				OpenGLColor curr_color = OpenGLColor(curr_rgb[0], curr_rgb[1], curr_rgb[2], 1.);
				multi_color.push_back(curr_color);
			}

			sas.clear();
			std::string sa_name = output_dir + "/" + std::to_string(frame) + "/" + sa_file_name;
			if (File::File_Exists(sa_name)) {
				BinaryDataIO::Read_Scalar_Array<real>(sa_name, sas);
				for (int i = 0; i < sas.size(); i++) {
					sas[i] = 200 * sqrt(sas[i]/3.14159);
				}
			}
			else {
				int pn = (int)point_vertices.size() / 8;
				sas.resize(pn);
				AuxFunc::Fill(sas, 1.);
			}
		}
		if (scale != last_scale) { Rescale(); }
		Set_Data_Refreshed();
	}

	void Update_Model_Matrix_Helper(const float* pos,const float* dir,const float radius,glm::mat4& model_matrix) const
	{
		model_matrix=glm::translate(glm::mat4(1.f),glm::vec3(pos[0],pos[1],pos[2]));
		float rad;Vector3f axis;AuxFunc::Angle_And_Axis_Between(Vector3f(0.f,0.f,1.f),Vector3f(dir[0],dir[1],dir[2]),rad,axis);
		model_matrix=glm::rotate(model_matrix,rad,glm::vec3(axis[0],axis[1],axis[2]));
		model_matrix=glm::scale(model_matrix,glm::vec3(radius,radius,radius));
	}

	virtual void Display() const
    {
		if (!visible) return;
		using namespace OpenGLUbos;

		switch(shading_mode){
		case ShadingMode::Multicolor: {
			std::shared_ptr<OpenGLShaderProgram> shader = shader_programs[0];
			GLuint stride_size = 4;
			shader->Begin();
			glPushAttrib(GL_LINE_BIT);
			Bind_Uniform_Block_To_Ubo(shader, "camera");
			glLineWidth(line_width);
			//shader->Set_Uniform_Vec4f("color", color.rgba);
			int pn = (int)point_vertices.size() / 8;
			for (int i = 0; i < pn; i++) {
				glm::mat4 model_p; Update_Model_Matrix_Helper(&point_vertices[i * 8], &point_vertices[i * 8 + 4], radius * scale * sas[i], model_p);
				shader->Set_Uniform_Matrix4f("model", glm::value_ptr(model_p));
				OpenGLColor color_code;
				color_code = multi_color[i];
				shader->Set_Uniform_Vec4f("color", color_code.rgba);
				glBindVertexArray(vao);
				glDrawArrays(GL_POLYGON, 0, vtx_size / stride_size);
			}
			glPopAttrib();
			shader->End();
		}break;
		case ShadingMode::Lighting:{
			std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[1];
			GLuint stride_size=8;

			shader->Begin();
			Bind_Uniform_Block_To_Ubo(shader,"camera");
			Bind_Uniform_Block_To_Ubo(shader,"lights");
			int pn=(int)point_vertices.size()/8;
			for(int i=0;i<pn;i++){
				glm::mat4 model_p;Update_Model_Matrix_Helper(&point_vertices[i*8],&point_vertices[i*8+4], radius * scale,model_p);
				shader->Set_Uniform_Matrix4f("model",glm::value_ptr(model_p));
				shader->Set_Uniform_Mat(&mat);
				shader->Set_Uniform_Vec4f("color", color.rgba);
				glBindVertexArray(vao);
				switch(polygon_mode){
				case PolygonMode::Fill:
					glDrawArrays(GL_POLYGON,0,vtx_size/stride_size);break;
				case PolygonMode::Wireframe:
					glDrawArrays(GL_LINE_LOOP,0,vtx_size/stride_size);break;}}
			shader->End();
		}break;
		case ShadingMode::ColorEncoding:{
			std::shared_ptr<OpenGLShaderProgram> shader = shader_programs[0];
			GLuint stride_size = 4;
			shader->Begin();
			glPushAttrib(GL_LINE_BIT);
			Bind_Uniform_Block_To_Ubo(shader, "camera");
			glLineWidth(line_width);
			//shader->Set_Uniform_Vec4f("color", color.rgba);
			int pn = (int)point_vertices.size() / 8;
			for (int i = 0; i < pn; i++) {
				glm::mat4 model_p; Update_Model_Matrix_Helper(&point_vertices[i * 8], &point_vertices[i * 8 + 4], radius * scale, model_p);
				shader->Set_Uniform_Matrix4f("model", glm::value_ptr(model_p));
				OpenGLColor color_code;
				if (i == color_code_selection) color_code = OpenGLColor::White();
				else color_code = ID_To_Color(i + color_code_offset);
				shader->Set_Uniform_Vec4f("color", color_code.rgba);
				glBindVertexArray(vao);
				glDrawArrays(GL_POLYGON, 0, vtx_size / stride_size);}
			glPopAttrib();
			shader->End();
		}break;
		default:{
			std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
			GLuint stride_size=4;

			shader->Begin();
			glPushAttrib(GL_LINE_BIT);
			Bind_Uniform_Block_To_Ubo(shader,"camera");
			glLineWidth(line_width);
			shader->Set_Uniform_Vec4f("color",color.rgba);	

			int pn=(int)point_vertices.size()/8;
			for (int i = 0; i < pn; i++) {
				glm::mat4 model_p; Update_Model_Matrix_Helper(&point_vertices[i * 8], &point_vertices[i * 8 + 4], radius * scale, model_p);
				shader->Set_Uniform_Matrix4f("model", glm::value_ptr(model_p));
				glBindVertexArray(vao);
				switch (polygon_mode) {
				case PolygonMode::Fill:
					glDrawArrays(GL_POLYGON, 0, vtx_size / stride_size); break;
				case PolygonMode::Wireframe:
					glDrawArrays(GL_LINE_LOOP, 0, vtx_size / stride_size); break;
				}
			}

			glPopAttrib();
			shader->End();
		}break;
		};
    }

	void Rescale()
	{
		//if (scale != last_scale) { std::cout << "new scale: " << scale << std::endl; }
		//float step = (float)two_pi / (float)n;
		//for (int i = 0; i < n; i++) {
		//	//OpenGL_Vertex4(Vector3(scale * cos((float)i * step), scale * sin((float)i * step), scale * (float)0), circle_vertices);	////position, 4 floats
		//	//OpenGL_Vertex4(Vector3(0., 0., 1.), circle_vertices);
		//	circle_vertices[8 * i + 0] = scale * cos((float)i * step);
		//	circle_vertices[8 * i + 1] = scale * sin((float)i * step);
		//}											////normal, 4 floats
		//vtx_size = circle_vertices.size();
		//Base::Set_OpenGL_Vertices(circle_vertices, vtx_size);
		//Set_OpenGL_Vertex_Attribute(0, 4, 8, 0);	////position
		//Set_OpenGL_Vertex_Attribute(1, 4, 8, 4);	////normal

		//last_scale = scale;
	}

	virtual int Color_Code_Size(void) {
		return (int)point_vertices.size() / 8;
	}

	virtual void Output_Debug_Info(std::ostream& out) {
		Base::Output_Debug_Info(out);
		out << ", size: " << Color_Code_Size();
	}

	void Toggle_Increase_Feature() {
		color_multiplier *= 1.01;
		std::cout << "increase color_multiplier to: " << color_multiplier << std::endl;
	}

	void Toggle_Decrease_Feature() {
		color_multiplier *= 0.99;
		std::cout << "decrease color_multiplier to: " << color_multiplier << std::endl;
	}

	void Set_Color_Multiplier(real cm) {
		color_multiplier = cm;
	}

	void Set_Radius(real multiplier)
	{
		radius = multiplier;
		//last_scale = scale;
	}
};
#endif