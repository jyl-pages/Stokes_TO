//////////////////////////////////////////////////////////////////////////
// Opengl particles
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLParticles_h__
#define __OpenGLParticles_h__
#include "Common.h"
#include "Particles.h"
#include "File.h"
#include "OpenGLObject.h"
#include "OpenGLVectors.h"
#include "OpenGLPoints.h"
#include "OpenGLMarkerObjects.h"

//////////////////////////////////////////////////////////////////////////
////The T_PARTICLE type must have attributes X, F, C

template<class T_PARTICLE=Particles<3> >
class OpenGLParticles : public OpenGLObject
{typedef OpenGLObject Base;
public:
    T_PARTICLE particles;
	
	OpenGLPoints opengl_points;
	Array<OpenGLVectors> opengl_vector_fields;
	OpenGLVectors* opengl_vector_velocity=nullptr;
	OpenGLVectors* opengl_vector_force=nullptr;
	
	ParticleMode particle_mode=ParticleMode::Dot;
	OpenGLCircle opengl_circle;
	OpenGLSphere opengl_sphere;
	bool use_particle_radius=false;

	virtual void Set_Particle_Mode(const ParticleMode _mode)
	{
		particle_mode=_mode;
		switch(particle_mode){
		case ParticleMode::Circle: opengl_circle.Initialize();break;
		case ParticleMode::Sphere: opengl_sphere.Initialize();break;}
	}

	virtual void Set_Particle_Size(const real p_size)
	{
		switch(particle_mode){
		case ParticleMode::Circle: {opengl_circle.radius=p_size;opengl_circle.Initialize();opengl_circle.Set_Data_Refreshed();opengl_circle.Update_Data_To_Render();}break;
		case ParticleMode::Sphere: {opengl_sphere.radius=p_size;opengl_sphere.Initialize();opengl_sphere.Set_Data_Refreshed();opengl_sphere.Update_Data_To_Render();}break;}
	}

	OpenGLParticles(){color=OpenGLColor::Red();name="particles";}

	virtual void Initialize()
	{
		Base::Initialize();
		opengl_points.Initialize();
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("particle"));
		T_PARTICLE* particles=nullptr;Initialize_Vector_Fields_Helper(particles);
	
		switch(particle_mode){
		case ParticleMode::Circle:
			opengl_circle.Initialize();break;
		case ParticleMode::Sphere:
			opengl_sphere.Initialize();break;}
	}

	void Set_Data_Pointers()
	{
		Set_Data_Pointers_Helper(particles);
	}

	virtual void Update_Data_To_Render()
	{
		if(!Update_Data_To_Render_Pre())return;
		Set_Data_Pointers();
		opengl_points.Update_Data_To_Render();
		for(auto& vf:opengl_vector_fields){vf.Update_Data_To_Render();}
		
		Update_Data_To_Render_Post();
	}

	virtual void Refresh(const int frame)
	{
		std::string file_name=output_dir+"/"+std::to_string(frame)+"/"+name;
		if(File::File_Exists(file_name)){
			File::Read_Binary_From_File(file_name,particles);
			Set_Data_Pointers();
			Set_Data_Refreshed();}
	}

	virtual void Display() const
    {
		if(!visible)return;
		switch(particle_mode){
		case ParticleMode::Dot:
			opengl_points.Display();
			break;
		case ParticleMode::Circle:
			if(use_particle_radius) opengl_circle.Display_Multiple_Instances(particles.XRef(),particles.CRef());
			else opengl_circle.Display_Multiple_Instances(particles.XRef());
			break;
		case ParticleMode::Sphere:
			opengl_sphere.Display_Multiple_Instances(particles.XRef());
			break;}

		for(auto& vf:opengl_vector_fields){vf.Display();}
    }

	virtual void Set_Color(const OpenGLColor& c)
	{
		color=c;opengl_points.color=c;opengl_circle.color=c;opengl_sphere.color=c;
	}

	virtual void Set_Shading_Mode(const ShadingMode _shading_mode)
	{shading_mode=_shading_mode;opengl_points.shading_mode=_shading_mode;}
	virtual void Set_Point_Size(const GLfloat point_size)
	{opengl_points.point_size=point_size;}
	virtual void Set_Point_Size(const Array<GLfloat>& point_size)
	{opengl_points.varying_point_size=point_size;opengl_points.use_varying_point_size=true;}
	virtual void Set_Line_Width(const GLfloat _line_width)
	{line_width=opengl_points.line_width=opengl_circle.line_width=opengl_sphere.line_width=_line_width;}

	void Set_Visibility_For_Velocity_Field(const bool vis=true){if(opengl_vector_velocity)opengl_vector_velocity->visible=vis;}
	void Set_Visibility_For_Force_Field(const bool vis=true){if(opengl_vector_force){opengl_vector_force->visible=vis;}}

	virtual void Toggle_Draw_Velocity(){if(opengl_vector_velocity)opengl_vector_velocity->visible=!opengl_vector_velocity->visible;}
    Define_Function_Object(OpenGLParticles<T_PARTICLE>,Toggle_Draw_Velocity);

	virtual void Toggle_Draw_Force(){if(opengl_vector_force)opengl_vector_force->visible=!opengl_vector_force->visible;}
	Define_Function_Object(OpenGLParticles<T_PARTICLE>,Toggle_Draw_Force);

	virtual void Toggle_Decrease_Scale()
	{Base::Toggle_Decrease_Scale();
	if(opengl_vector_velocity)opengl_vector_velocity->Toggle_Decrease_Scale();
	if(opengl_vector_force)opengl_vector_force->Toggle_Decrease_Scale();}
    Define_Function_Object(OpenGLParticles<T_PARTICLE>,Toggle_Decrease_Scale);

	virtual void Toggle_Increase_Scale()
	{Base::Toggle_Increase_Scale();
	if(opengl_vector_velocity)opengl_vector_velocity->Toggle_Increase_Scale();
	if(opengl_vector_force)opengl_vector_force->Toggle_Increase_Scale();}
    Define_Function_Object(OpenGLParticles<T_PARTICLE>,Toggle_Increase_Scale);

protected:
	void Initialize_Vector_Fields_Helper(Points<3>* particles=nullptr){}
	void Initialize_Vector_Fields_Helper(Particles<3>* particles=nullptr)
	{opengl_vector_fields.resize(2);for(auto& vf:opengl_vector_fields){vf.Initialize();vf.visible=false;/*vector fields invisble by default*/}}

	void Set_Data_Pointers_Helper(Points<3>& particles)
	{opengl_points.Set_Data_Pointers(particles.X());}
	void Set_Data_Pointers_Helper(Particles<3>& particles)
	{
		if(!initialized)Initialize();
		opengl_points.Set_Data_Pointers(particles.X(),particles.C());
		opengl_vector_fields[0].Set_Data_Pointers(particles.X(),particles.V());
		opengl_vector_velocity=&opengl_vector_fields[0];
		opengl_vector_fields[1].Set_Data_Pointers(particles.X(),particles.F());
		opengl_vector_force=&opengl_vector_fields[1];
	}
};
#endif