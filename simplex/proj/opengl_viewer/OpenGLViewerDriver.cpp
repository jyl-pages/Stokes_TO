//////////////////////////////////////////////////////////////////////////
// Opengl viewer driver
// Copyright (c) (2018-),Bo Zhu
// This file is part of SimpleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "OpenGLViewerDriver.h"
#include "OpenGLViewer.h"
#include "OpenGLGrid.h"
#include "OpenGLGridField.h"
#include "OpenGLGridVectors.h"
#include "OpenGLGridTensors.h"
#include "OpenGLGridHeightField.h"
#include "OpenGLMesh.h"
#include "OpenGLMeshField.h"
#include "OpenGLMeshVectors.h"
#include "OpenGLMeshTensors.h"
#include "OpenGLParticles.h"
#include "OpenGLTrackerPoints.h"
#include "OpenGLTrackerVectors.h"
#include "OpenGLTrackerCircles.h"
#include "OpenGLScreenObjects.h"
#include "OrientedParticle.h"
#include "OpenGLAabb.h"
#include "OpenGLWindow.h"

#include <iostream>

//#include "DataIO.h"

//// pybind
//#include <pybind11/pybind11.h>
//#include <pybind11/embed.h>
//#include <pybind11/numpy.h>
//#include <pybind11/eigen.h>
//#include <pybind11/stl.h>    // for myPyObject.cast<std::vector<T>>()
//namespace py=pybind11;

//////////////////////////////////////////////////////////////////////////
////fluid viewers

void OpenGLViewerFluidEuler::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::Gray(.4f),OpenGLColor::Gray(.2f));
}

void OpenGLViewerFluidEuler::Initialize_Data()
{
	////grid values
	auto grid=Add_Object<OpenGLGrid>("grid");Set_Visibility(grid,'g',true);
	//Set_Subgrid_Slice(grid,2,10);

	auto vec=Add_Grid_Vector_Field("grid",{{"velocity","1",ColorType::Jet,StoreType::Cell}});Set_Visibility(vec,'v',true);
	Set_Auto_Length_Norm(vec,1.);
	//Set_Subgrid_Slice(vec,2,10);


	auto phi=Add_Grid_Scalar_Field("grid",{{"phi","",ColorType::Jet,StoreType::Cell}});Set_Visibility(phi,'l',false);
	//auto phi_hf=Add_Grid_Height_Field("grid",{{"phi","",ColorType::Jet,StoreType::Node}});Set_Visibility(phi,'L',false);
	auto vor=Add_Grid_Scalar_Field("grid",{{"vorticity","",ColorType::Jet,StoreType::Cell}});Set_Visibility(vor,'o',false);
	auto type=Add_Grid_Scalar_Field("grid",{{"fluid_type","",ColorType::Jet,StoreType::Cell}});Set_Visibility(type,'T',false);
	auto f=Add_Grid_Vector_Field("grid",{{"force","1",ColorType::Jet,StoreType::Cell}});Set_Visibility(f,'F',false);
	auto a=Add_Grid_Vector_Field("grid",{{"acc","1",ColorType::Jet,StoreType::Cell}});Set_Visibility(a,'a',false);

	////particle
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("particles");Set_Visibility(p,'f',true);
	Set_Color(p,OpenGLColor::Green());Set_Point_Size(p,2.f);
	Set_Particle_Mode(p,ParticleMode::Dot);/*Set_Particle_Size(p,.1f);*/}

	////interface
	auto c2=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(c2,'s',true);Set_Color(c2,OpenGLColor::Red());Set_Line_Width(c2,4.f);
	auto c3=Add_Object<OpenGLTriangleMesh>("triangle_mesh");Set_Polygon_Mode(c3,PolygonMode::SurfOnly);Set_Color(c3,OpenGLColor::Red());Set_Visibility(c3,'t',true);
	auto c4=Add_Mesh_Scalar_Field(c3,{{"triangle_mesh_node_color","",ColorType::Jet,StoreType::Node}});Set_Visibility(c4,'c',true);
	auto c5=Add_Object<OpenGLTetrahedronMesh>("tetrahedron_mesh");Set_Polygon_Mode(c5,PolygonMode::Wireframe);Set_Color(c5,OpenGLColor::Red());Set_Visibility(c5,'T',true);

	////boundary
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("psi_D");Set_Visibility(p,'d',true);Set_Color(p,OpenGLColor::Red());}
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("psi_N");Set_Visibility(p,'n',true);Set_Color(p,OpenGLColor::Blue());}
}

void OpenGLViewerFluidEulerHighRes::Initialize_Data()
{
	auto grid=Add_Object<OpenGLGrid>("grid");Set_Visibility(grid,'g',true);
	auto type=Add_Grid_Scalar_Field("grid",{{"fluid_type","",ColorType::Jet,StoreType::Cell}});Set_Visibility(type,'T',false);
	auto vec=Add_Grid_Vector_Field("grid",{{"velocity","1",ColorType::Jet,StoreType::Cell}});Set_Visibility(vec,'v',true);

	{auto p=Add_Object<OpenGLTrackerPoints>("tracker_points");Set_Visibility(p,'f',true);
	Set_Color(p,OpenGLColor::Green());Set_Point_Size(p,1.f);}
}

void OpenGLViewerFluidLagrangian::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::Gray(.4f),OpenGLColor::Gray(.2f));
}

void OpenGLViewerFluidLagrangian::Initialize_Data()
{
	////grid
	auto grid=Add_Object<OpenGLGrid>("grid");Set_Visibility(grid,'g',true);
	{auto vec=Add_Grid_Vector_Field("grid",{ {"velocity","",ColorType::Jet,StoreType::Cell}});Set_Visibility(vec,'v',true);Set_Scale(vec,default_scale);}
	{auto vec=Add_Grid_Vector_Field("grid",{ {"velocity2","",ColorType::Jet,StoreType::Cell}});Set_Visibility(vec,'z',true);Set_Scale(vec,default_scale);}
	{auto vec=Add_Grid_Vector_Field("grid",{ {"impulse","",ColorType::Jet,StoreType::Cell}});Set_Visibility(vec,'a',true);Set_Scale(vec,default_scale);}
	auto type=Add_Grid_Scalar_Field("grid",{ {"fluid_type","",ColorType::Jet,StoreType::Cell}});Set_Visibility(type,'T',false);
	{auto psi=Add_Grid_Vector_Field("grid",{ {"psi","",ColorType::Jet,StoreType::Cell}});Set_Visibility(psi,'S',true);Set_Scale(psi,default_scale);}

	auto phi=Add_Grid_Scalar_Field("grid",{ {"phi","",ColorType::Jet,StoreType::Cell}});Set_Visibility(phi,'l',false);

	////particle
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("particles");Set_Visibility(p,'f',true);
	Set_Color(p,OpenGLColor::Red());Set_Point_Size(p,8.f);Set_Particle_Velocity_Visibility(p,'v');Set_Particle_Force_Visibility(p,'F');}

	{auto p=Add_Object<OpenGLTrackerPoints>("tracker_points");Set_Visibility(p,'T',true);Set_Color(p,OpenGLColor::Yellow());Set_Point_Size(p,4.f);}

	{auto p=Add_Object<OpenGLTrackerCircles>("tracker_circles");Set_Visibility(p,'f',false);
	Set_Color(p,OpenGLColor::Yellow());Set_Point_Size(p,1.f);Set_Shading_Mode(p,ShadingMode::Lighting);}

	{auto p=Add_Object<OpenGLTrackerVectors>("point_force");Set_Visibility(p,'F',true);Set_Color(p,OpenGLColor::Green());Set_Scale(p,default_scale);}

	{auto p=Add_Object<OpenGLTrackerVectors>("point_velocity");Set_Visibility(p,'V',true);Set_Color(p,OpenGLColor::Blue());Set_Scale(p,default_scale);}

	{auto p=Add_Object<OpenGLTrackerVectors>("point_impulse");Set_Visibility(p,'I',true);Set_Color(p,OpenGLColor::Red());Set_Scale(p,default_scale);}

	{auto p=Add_Object<OpenGLTrackerVectors>("point_gamma");Set_Visibility(p,'G',false);Set_Color(p,OpenGLColor::Yellow());Set_Scale(p,default_scale);}

	{auto p=Add_Object<OpenGLTrackerVectors>("point_height");Set_Visibility(p,'H',false);Set_Color(p,OpenGLColor::Red());Set_Scale(p,default_scale);}

	////mesh
	auto c2=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(c2,'s',true);Set_Color(c2,OpenGLColor::Red());Set_Line_Width(c2,4.f);	////for levelset contour

	//auto c2=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(c2,'s',false);Set_Color(c2,OpenGLColor::Green());Set_Line_Width(c2,1.f);	////for local frame
	//if(!use_2d_display){
	//	Array<float> multi_color={1.f,0.f,0.f,1.f,1.f,0.f,0.f,1.f,0.f,1.f,0.f,1.f,0.f,1.f,0.f,1.f,0.f,0.f,1.f,1.f,0.f,0.f,1.f,1.f};	////RGB for three axes
	//	Set_Multi_Color(c2,multi_color);}
	//else{
	//	Array<float> multi_color={1.f,0.f,0.f,1.f,1.f,0.f,0.f,1.f,0.f,1.f,0.f,1.f,0.f,1.f,0.f,1.f};	////RGB for three axes
	//	Set_Multi_Color(c2,multi_color);}

	//auto c4=Add_Object<OpenGLSegmentMesh>("segment_mesh2");Set_Visibility(c4,'F',true);Set_Color(c4,OpenGLColor::Green());Set_Line_Width(c4,1.f);

	{auto p=Add_Object<OpenGLTrackerVectors>("tracker_vectors");Set_Visibility(p,'F',true);Set_Color(p,OpenGLColor::Green());Set_Point_Size(p,1.f);}

	auto c3=Add_Object<OpenGLTriangleMesh>("triangle_mesh");Set_Polygon_Mode(c3,PolygonMode::Wireframe);Set_Shading_Mode(c3,ShadingMode::None);
	Set_Color(c3,OpenGLColor::Green(.5f));Set_Visibility(c3,'t',true);Set_Line_Width(c3,4.f);

	////boundary
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("psi_D");Set_Visibility(p,'d',false);Set_Color(p,OpenGLColor::Red());}
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("psi_N");Set_Visibility(p,'n',false);Set_Color(p,OpenGLColor::Blue());}

	////lighting
	//auto dir_light=OpenGLUbos::Add_Directional_Light(glm::vec3(-1.f,-.1f,-.2f));
	////dir_light->dif=glm::vec4(1.,0.,0.,1.);
	//OpenGLUbos::Set_Ambient(glm::vec4(.1f,.1f,.1f,1.f));
	//OpenGLUbos::Update_Lights_Ubo();
}


void OpenGLViewerFluidSPHBubble::Initialize()
{
	Base::Initialize();
	//Set_Background_Color(OpenGLColor::Black(),OpenGLColor::Black());
	Set_Background_Color(OpenGLColor::White(), OpenGLColor::White());
}

void OpenGLViewerFluidSPHBubble::Initialize_Data()
{
	////grid
	auto grid=Add_Object<OpenGLGrid>("grid");Set_Visibility(grid,'g',true);
	{auto vec=Add_Grid_Vector_Field("grid",{ {"velocity","",ColorType::Jet,StoreType::Cell}});Set_Visibility(vec,'v',false);Set_Scale(vec,default_scale);}
	{auto vec=Add_Grid_Vector_Field("grid",{ {"velocity2","",ColorType::Jet,StoreType::Cell}});Set_Visibility(vec,'z',false);Set_Scale(vec,default_scale);}	
	{auto vec=Add_Grid_Vector_Field("grid",{ {"impulse","",ColorType::Jet,StoreType::Cell}});Set_Visibility(vec,'a',false);Set_Scale(vec,default_scale);}	
	auto type=Add_Grid_Scalar_Field("grid",{{"fluid_type","",ColorType::Jet,StoreType::Cell}});Set_Visibility(type,'T',false);
	
	auto phi=Add_Grid_Scalar_Field("grid",{{"phi","",ColorType::Jet,StoreType::Cell}});Set_Visibility(phi,'l',false);

	////particle
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("particles");Set_Visibility(p,'f',false);
	Set_Color(p,OpenGLColor::Red());Set_Point_Size(p,8.f);Set_Particle_Velocity_Visibility(p,'v');Set_Particle_Force_Visibility(p,'F');}

	{auto p=Add_Object<OpenGLTrackerPoints>("tracker_points");Set_Visibility(p,'T',true);Set_Color(p,OpenGLColor::Orange());Set_Point_Size(p,4.f);}
	{auto p = Add_Object<OpenGLTrackerPoints>("tracker_points_3d"); Set_Visibility(p, 'T', true); Set_Color(p, OpenGLColor::Green()); Set_Point_Size(p, 2.f); }

	{auto p=Add_Object<OpenGLTrackerCircles>("tracker_circles");Set_Visibility(p,'f',true);
	p->Set_Height_File("h_bin");
	Set_Shading_Mode(p,ShadingMode::Multicolor);
	}

	{auto p=Add_Object<OpenGLTrackerVectors>("point_force");Set_Visibility(p,'F',false);Set_Color(p,OpenGLColor::Green());Set_Scale(p,default_scale);}
	
	{auto p=Add_Object<OpenGLTrackerVectors>("point_velocity");Set_Visibility(p,'V',false);Set_Color(p,OpenGLColor::Blue());Set_Scale(p,default_scale);}

	{auto p=Add_Object<OpenGLTrackerVectors>("point_impulse");Set_Visibility(p,'I',false);Set_Color(p,OpenGLColor::Red());Set_Scale(p,default_scale);}

	{auto p=Add_Object<OpenGLTrackerVectors>("point_gamma");Set_Visibility(p,'G',false);Set_Color(p,OpenGLColor::Yellow());Set_Scale(p,default_scale);}

	{auto p=Add_Object<OpenGLTrackerVectors>("point_height");Set_Visibility(p,'H',false);Set_Color(p,OpenGLColor::Red());Set_Scale(p,default_scale);}

	////mesh
	auto c2=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(c2,'s',true);Set_Color(c2,OpenGLColor::Red());Set_Line_Width(c2,2.f);	////for levelset contour

	//auto c2=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(c2,'s',false);Set_Color(c2,OpenGLColor::Green());Set_Line_Width(c2,1.f);	////for local frame
	//if(!use_2d_display){
	//	Array<float> multi_color={1.f,0.f,0.f,1.f,1.f,0.f,0.f,1.f,0.f,1.f,0.f,1.f,0.f,1.f,0.f,1.f,0.f,0.f,1.f,1.f,0.f,0.f,1.f,1.f};	////RGB for three axes
	//	Set_Multi_Color(c2,multi_color);}
	//else{
	//	Array<float> multi_color={1.f,0.f,0.f,1.f,1.f,0.f,0.f,1.f,0.f,1.f,0.f,1.f,0.f,1.f,0.f,1.f};	////RGB for three axes
	//	Set_Multi_Color(c2,multi_color);}

	//auto c4=Add_Object<OpenGLSegmentMesh>("segment_mesh2");Set_Visibility(c4,'F',true);Set_Color(c4,OpenGLColor::Green());Set_Line_Width(c4,1.f);

	{auto p=Add_Object<OpenGLTrackerVectors>("tracker_vectors");Set_Visibility(p,'F',false);Set_Color(p,OpenGLColor::Green());Set_Point_Size(p,1.f);}

	auto c3=Add_Object<OpenGLTriangleMesh>("triangle_mesh");Set_Polygon_Mode(c3,PolygonMode::Wireframe);Set_Shading_Mode(c3,ShadingMode::None);
	Set_Color(c3,OpenGLColor::Green(.5f));Set_Visibility(c3,'t',false);Set_Line_Width(c3,4.f);
	
	////boundary
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("psi_D");Set_Visibility(p,'d',false);Set_Color(p,OpenGLColor::Red());}
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("psi_N");Set_Visibility(p,'n',false);Set_Color(p,OpenGLColor::Blue());}

	////lighting
	//auto dir_light=OpenGLUbos::Add_Directional_Light(glm::vec3(-1.f,-.1f,-.2f));
	////dir_light->dif=glm::vec4(1.,0.,0.,1.);
	//OpenGLUbos::Set_Ambient(glm::vec4(.1f,.1f,.1f,1.f));
	//OpenGLUbos::Update_Lights_Ubo();
}

void OpenGLViewerParticleALEFilm::Initialize()
{
	Base::Initialize();
	//Set_Background_Color(OpenGLColor::Black(),OpenGLColor::Black());
	Set_Background_Color(OpenGLColor::Black(), OpenGLColor::Black());
}

void OpenGLViewerParticleALEFilm::Initialize_Data()
{
	{auto p = Add_Object<OpenGLTrackerPoints>("boundary_points");Set_Visibility(p, 'b', true);Set_Color(p, OpenGLColor::White());Set_Point_Size(p, 4.f);}
	{auto p = Add_Object<OpenGLTrackerPoints>("e_tracker_points");Set_Visibility(p, 'T', true);Set_Color(p, OpenGLColor::Orange());Set_Point_Size(p, 4.f);}
	{auto p = Add_Object<OpenGLTrackerCircles>("e_tracker_circles");Set_Visibility(p, 'C', false);
	Set_Shading_Mode(p, ShadingMode::Multicolor);
	if (p)p->Set_Height_File("e_h_bin");
	if (p)p->Set_Radius(0.03);
	if (p)p->Set_Color_Multiplier(1e9);
	}
	{auto p = Add_Object<OpenGLTrackerVectors>("e_point_force");Set_Visibility(p, 'F', false);Set_Color(p, OpenGLColor::Green());Set_Scale(p, 20 * default_scale);}
	{auto p = Add_Object<OpenGLTrackerVectors>("e_point_velocity");Set_Visibility(p, 'V', false);Set_Color(p, OpenGLColor::Blue());Set_Scale(p, default_scale);}
	{auto p = Add_Object<OpenGLTrackerVectors>("e_point_height");Set_Visibility(p, 'H', false);Set_Color(p, OpenGLColor::Red());Set_Scale(p, default_scale);}
	{auto p = Add_Object<OpenGLTrackerVectors>("e_point_nden");Set_Visibility(p, 'J', false);Set_Color(p, OpenGLColor::Yellow());Set_Scale(p, default_scale);}
	{auto p = Add_Object<OpenGLTrackerPoints>("e_point_temperature"); Set_Visibility(p, 'e', false); Set_Color(p, OpenGLColor::Red()); }
	{auto p = Add_Object<OpenGLTrackerVectors>("e_point_pressure");Set_Visibility(p, 'P', false);Set_Color(p, OpenGLColor::Yellow());Set_Scale(p, default_scale);}
	{auto p = Add_Object<OpenGLTrackerVectors>("e_point_h_normal");Set_Visibility(p, 'N', false);Set_Color(p, OpenGLColor::Gray());Set_Scale(p, default_scale);}
	{auto p = Add_Object<OpenGLTrackerVectors>("e_point_3d_force");Set_Visibility(p, 'I', false);Set_Color(p, OpenGLColor::Orange());Set_Scale(p, default_scale);}
	{auto p = Add_Object<OpenGLTrackerVectors>("e_point_h_curvature");Set_Visibility(p, 'K', false);Set_Color(p, OpenGLColor::Blue());Set_Scale(p, default_scale);}

	{auto p = Add_Object<OpenGLTrackerVectors>("e_point_IPN");Set_Visibility(p, 'Z', false);Set_Color(p, OpenGLColor::Green());Set_Scale(p, default_scale);}
	{auto p = Add_Object<OpenGLTrackerVectors>("e_point_IPK");Set_Visibility(p, 'X', false);Set_Color(p, OpenGLColor::Blue());Set_Scale(p, default_scale);}
	{auto p = Add_Object<OpenGLTrackerVectors>("e_point_IPF");Set_Visibility(p, 'x', false);Set_Color(p, OpenGLColor::White());Set_Scale(p, default_scale);}
	{auto p = Add_Object<OpenGLTrackerVectors>("e_point_CF");Set_Visibility(p, 'z', false);Set_Color(p, OpenGLColor::Orange());Set_Scale(p, default_scale);}

	{auto p = Add_Object<OpenGLTrackerVectors>("e_point_target_height");Set_Visibility(p, 'G', false);Set_Color(p, OpenGLColor::Yellow());Set_Scale(p, default_scale);}
	{auto p = Add_Object<OpenGLSegmentMesh>("e_segment_mesh");Set_Visibility(p, 'S', false);Set_Color(p, OpenGLColor::Red());Set_Line_Width(p, 2.f);}

	{auto p = Add_Object<OpenGLTrackerPoints>("l_tracker_points");Set_Visibility(p, 't', false);Set_Color(p, OpenGLColor::Green());Set_Point_Size(p, 4.f);}
	//{auto p = Add_Object<OpenGLTrackerCircles>("l_tracker_circles");Set_Visibility(p, 'c', true);
	//Set_Color(p, OpenGLColor::Yellow());
	//Set_Shading_Mode(p, ShadingMode::Multicolor);
	//}
	{auto p = Add_Object<OpenGLTrackerVectors>("l_point_force");Set_Visibility(p, 'f', false);Set_Color(p, OpenGLColor::Green());Set_Scale(p, default_scale);}
	{auto p = Add_Object<OpenGLTrackerVectors>("l_point_velocity");Set_Visibility(p, 'v', false);Set_Color(p, OpenGLColor::Blue());Set_Scale(p, default_scale);}
	{auto p = Add_Object<OpenGLTrackerVectors>("l_point_height");Set_Visibility(p, 'h', false);Set_Color(p, OpenGLColor::Red());Set_Scale(p, default_scale);}
	//{auto p = Add_Object<OpenGLSegmentMesh>("l_segment_mesh");Set_Visibility(p, 's', true);Set_Color(p, OpenGLColor::Blue());Set_Line_Width(p, 2.f);}
	{auto p = Add_Object<OpenGLTrackerCircles>("l_tracker_circles");Set_Visibility(p, 'c', true);
	Set_Shading_Mode(p, ShadingMode::Multicolor);
	if (p)p->Set_Height_File("l_h_bin");
	if (p)p->Set_SA_File("l_sa_bin");
	if (p)p->Set_Radius(0.01);
	if (p)p->Set_Color_Multiplier(1e9);
	}

	{auto p = Add_Object<OpenGLTrackerPoints>("v_tracker_points");Set_Visibility(p, '3', true);Set_Color(p, OpenGLColor::Blue());Set_Point_Size(p, 4.f);}

	//air solver part
	{auto grid = Add_Object<OpenGLGrid>("grid"); Set_Visibility(grid, 'g', true); }
	{auto vec = Add_Grid_Vector_Field("grid", { {"velocity","",ColorType::Jet,StoreType::Cell} }); Set_Visibility(vec, '1', true); Set_Scale(vec, default_scale); }
	{auto p = Add_Grid_Scalar_Field("grid", { {"temperature","",ColorType::Jet,StoreType::Cell} }); Set_Visibility(p, '2', false); }
	{auto p = Add_Grid_Scalar_Field("grid", { {"pressure","",ColorType::Jet,StoreType::Cell} }); Set_Visibility(p, '3', false); }
	{auto p = Add_Grid_Scalar_Field("grid", { {"phi","",ColorType::Jet,StoreType::Cell} }); Set_Visibility(p, '4', false); }
	{auto p = Add_Object<OpenGLParticles<Particles<3> > >("psi_D"); Set_Visibility(p, 'd', true); Set_Color(p, OpenGLColor::Red()); }
	{auto p = Add_Object<OpenGLParticles<Particles<3> > >("psi_N"); Set_Visibility(p, 'n', true); Set_Color(p, OpenGLColor::Blue()); }
}

void OpenGLViewerOrientedParticle::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::Gray(.4f),OpenGLColor::Gray(.2f));
}

void OpenGLViewerOrientedParticle::Initialize_Data()
{
	auto p=new OrientedParticle();
	p->Initialize();
	p->name="particles";
	p->pos.Initialize("x_bin",true,OpenGLColor::Yellow());
	p->normal.Initialize("norm_bin",false,OpenGLColor::Blue());
	p->Initialize_Data(output_dir,0);
	Add_OpenGL_Object(p);
	std::cout << "initialize data done with " << opengl_window->object_list.size() << std::endl;
}

void OpenGLViewerPoisson::Initialize_Data()
{
	auto c2=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(c2,'s',true);Set_Color(c2,OpenGLColor::Red());Set_Line_Width(c2,1.f);
	auto grid=Add_Object<OpenGLGrid>("grid");Set_Visibility(grid,'g',true);
	auto x=Add_Grid_Scalar_Field("grid",{{"x","1",ColorType::Jet,StoreType::Cell}});Set_Visibility(x,'x',true);
	//auto y = Add_Grid_Volume("grid", { {"x","",ColorType::Jet,StoreType::Cell} }); Set_Visibility(y, 'y', false);
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("psi_D");Set_Visibility(p,'d',true);Set_Color(p,OpenGLColor::Red());}
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("psi_N");Set_Visibility(p,'n',true);Set_Color(p,OpenGLColor::Blue());}
}

void OpenGLViewerLevelSet::Initialize_Data()
{
	auto grid=Add_Object<OpenGLGrid>("grid");Set_Visibility(grid,'g',true);
	auto phi=Add_Grid_Vector_Field("grid",{{"phi","",ColorType::Jet,StoreType::Cell}});Set_Visibility(phi,'p',false);
	auto curv=Add_Grid_Vector_Field("grid",{{"curvature","",ColorType::Jet,StoreType::Cell}});Set_Visibility(phi,'c',false);
	auto normal=Add_Grid_Vector_Field("grid",{{"normal","1",ColorType::Jet,StoreType::Cell}});Set_Visibility(normal,'n',true);

	auto c2=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(c2,'s',true);Set_Color(c2,OpenGLColor::Red());Set_Line_Width(c2,4.f);
	auto c3=Add_Object<OpenGLTriangleMesh>("triangle_mesh");Set_Polygon_Mode(c3,PolygonMode::Wireframe);Set_Visibility(c3,'t',true);
}

void OpenGLViewerFluidMesh::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::Gray(.4f),OpenGLColor::Gray(.2f));
}

void OpenGLViewerFluidMesh::Initialize_Data()
{
	auto p=Add_Object<OpenGLParticles<Particles<3> > >("particles");Set_Visibility(p,'x',true);Set_Color(p,OpenGLColor::Yellow());
	Set_Particle_Velocity_Visibility(p,'v');Set_Particle_Force_Visibility(p,'f');
	auto seg=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(seg,'s',true);Set_Color(seg,OpenGLColor::Red());
	auto tri=Add_Object<OpenGLTriangleMesh>("triangle_mesh");Set_Polygon_Mode(tri,PolygonMode::Wireframe);Set_Visibility(tri,'t',true);Set_Color(tri,OpenGLColor::Red());
	auto tet=Add_Object<OpenGLTetrahedronMesh>("tetrahedron_mesh");Set_Polygon_Mode(tet,PolygonMode::Wireframe);Set_Visibility(tet,'e',true);Set_Color(tet,OpenGLColor::Blue());
	
	if(tri!=nullptr){ 
		auto pressure=Add_Mesh_Scalar_Field(tri,{{"pressure","",ColorType::Jet,StoreType::Node}});Set_Visibility(pressure,'6',true);Set_Alpha(pressure,.5f);}
}

void OpenGLViewerMatchStick::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::Gray(.4f),OpenGLColor::Gray(.2f));
}

void OpenGLViewerMatchStick::Initialize_Data()
{
	auto grid=Add_Object<OpenGLGrid>("grid");Set_Visibility(grid,'g',true);
	{auto vec=Add_Grid_Vector_Field("grid",{ {"velocity","1",ColorType::Jet,StoreType::Cell},{"vorticity","2",ColorType::Jet,StoreType::Cell}});Set_Visibility(vec,'v',true);Set_Scale(vec,default_scale);}
	auto type=Add_Grid_Scalar_Field("grid",{{"fluid_type","",ColorType::Jet,StoreType::Cell}});Set_Visibility(type,'T',false);

	{auto p=Add_Object<OpenGLTrackerPoints>("points_red");Set_Visibility(p,'T',true);Set_Color(p,OpenGLColor::Red());Set_Point_Size(p,6.f);}
	{auto p=Add_Object<OpenGLTrackerPoints>("points_blue");Set_Visibility(p,'T',true);Set_Color(p,OpenGLColor::Blue());Set_Point_Size(p,6.f);}
	{auto p=Add_Object<OpenGLTrackerVectors>("sticks");Set_Visibility(p,'F',true);Set_Color(p,OpenGLColor::Green());Set_Scale(p,default_scale);}

	{auto phi=Add_Grid_Scalar_Field("grid",{{"phi","",ColorType::Jet,StoreType::Cell}});Set_Visibility(phi,'l',false);}

	{auto seg=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(seg,'s',true);Set_Color(seg,OpenGLColor::Red());}
	{auto tri=Add_Object<OpenGLTriangleMesh>("triangle_mesh");Set_Polygon_Mode(tri,PolygonMode::Wireframe);Set_Visibility(tri,'t',true);Set_Color(tri,OpenGLColor::Red());}
}

void OpenGLViewerFluidDEC::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::Gray(.4f), OpenGLColor::Gray(.2f));

}

void OpenGLViewerFluidDEC::Initialize_Data()
{
	////grid values
	auto grid = Add_Object<OpenGLGrid>("grid"); Set_Visibility(grid, 'g', true);
	auto tri = Add_Object<OpenGLTriangleMesh>("triangle_mesh"); Set_Polygon_Mode(tri, PolygonMode::Wireframe); Set_Visibility(tri, 't', true); Set_Color(tri, OpenGLColor::Red());
	
	if (grid != nullptr)
	{
		auto vec = Add_Grid_Vector_Field("grid", { {"velocity","1",ColorType::Jet,StoreType::Cell} }); Set_Visibility(vec, 'v', true);
		Set_Auto_Length_Norm(vec, 1.);

		auto phi = Add_Grid_Scalar_Field("grid", { {"phi","",ColorType::Jet,StoreType::Cell} }); Set_Visibility(phi, 'l', false);
		auto vor = Add_Grid_Scalar_Field("grid", { {"vorticity","",ColorType::Jet,StoreType::Cell} }); Set_Visibility(vor, 'o', false);
		auto type = Add_Grid_Scalar_Field("grid", { {"fluid_type","",ColorType::Jet,StoreType::Cell} }); Set_Visibility(type, 'T', false);
		auto f = Add_Grid_Vector_Field("grid", { {"force","1",ColorType::Jet,StoreType::Cell} }); Set_Visibility(f, 'F', false);
		auto a = Add_Grid_Vector_Field("grid", { {"acc","1",ColorType::Jet,StoreType::Cell} }); Set_Visibility(a, 'a', false);
		auto pressure = Add_Grid_Scalar_Field("grid", { {"pressure","1",ColorType::Jet,StoreType::Cell} }); Set_Visibility(pressure, 'P', false);

		auto x = Add_Grid_Scalar_Field("grid", { {"x","",ColorType::Den,StoreType::Cell} }); Set_Visibility(x, 'x', false); Set_Color(x, OpenGLColor::Red());
		auto y = Add_Grid_Volume("grid", { {"x","",ColorType::Den,StoreType::Cell} }); Set_Visibility(y, 'y', false);
	}

	if (tri != nullptr)
	{
		auto x = Add_Mesh_Scalar_Field(tri, { {"x","",ColorType::Jet,StoreType::Cell} }); Set_Visibility(x, 'x', true);
	}

}

void OpenGLViewerSWE::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::Gray(.4f), OpenGLColor::Gray(.2f));
}

void OpenGLViewerSWE::Initialize_Data()
{
	////grid values
	auto grid = Add_Object<OpenGLGrid>("grid"); Set_Visibility(grid, 'g', true);

	{auto vec = Add_Grid_Vector_Field("grid", { {"velocity","1",ColorType::Jet,StoreType::Cell} }); Set_Visibility(vec, 'v', true);
	Set_Auto_Length_Norm(vec, 1.);}

	{auto imp = Add_Grid_Vector_Field("grid", { {"impulse","1",ColorType::Jet,StoreType::Cell} }); Set_Visibility(imp, 'a', true);
	Set_Auto_Length_Norm(imp, 1.);}

	////particle
	{auto p = Add_Object<OpenGLParticles<Particles<3> > >("particles"); Set_Visibility(p, 'f', true);
	Set_Color(p, OpenGLColor::Green()); Set_Point_Size(p, 2.f);
	Set_Particle_Mode(p, ParticleMode::Dot);/*Set_Particle_Size(p,.1f);*/}

	////interface
	auto c2 = Add_Object<OpenGLSegmentMesh>("segment_mesh"); Set_Visibility(c2, 's', true); Set_Color(c2, OpenGLColor::Red()); Set_Line_Width(c2, 4.f);
	auto c3 = Add_Object<OpenGLTriangleMesh>("triangle_mesh"); Set_Polygon_Mode(c3, PolygonMode::SurfOnly); Set_Color(c3, OpenGLColor::Red()); Set_Visibility(c3, 't', true);
	auto c4 = Add_Mesh_Scalar_Field(c3, { {"triangle_mesh_node_color","",ColorType::Jet,StoreType::Node} }); Set_Visibility(c4, 'c', true);
	auto c5 = Add_Object<OpenGLTetrahedronMesh>("tetrahedron_mesh"); Set_Polygon_Mode(c5, PolygonMode::Wireframe); Set_Color(c5, OpenGLColor::Red()); Set_Visibility(c5, 'T', true);

	////boundary
	{auto p = Add_Object<OpenGLParticles<Particles<3> > >("psi_D"); Set_Visibility(p, 'd', true); Set_Color(p, OpenGLColor::Red()); }
	{auto p = Add_Object<OpenGLParticles<Particles<3> > >("psi_N"); Set_Visibility(p, 'n', true); Set_Color(p, OpenGLColor::Blue()); }
}

void OpenGLViewerMicroFluidic::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::Gray(.4f), OpenGLColor::Gray(.2f));
}

void OpenGLViewerMicroFluidic::Initialize_Data()
{
	////grid values
	auto grid = Add_Object<OpenGLGrid>("grid"); Set_Visibility(grid, 'g', true);
	auto type = Add_Grid_Scalar_Field("grid", { {"fluid_type","",ColorType::Jet,StoreType::Cell} }); Set_Visibility(type, 'T', false);
	auto vec = Add_Grid_Vector_Field("grid", { {"velocity","1",ColorType::Jet,StoreType::Cell} }); Set_Visibility(vec, 'v', true);
	Set_Auto_Length_Norm(vec, 1.);

	{auto p = Add_Object<OpenGLTrackerPoints>("tracker_points"); Set_Visibility(p, 'f', true);
	Set_Color(p, OpenGLColor::Green()); Set_Point_Size(p, 4.f); }
	{auto p = Add_Object<OpenGLTrackerPoints>("trapped_points"); Set_Visibility(p, 'F', true);
	Set_Color(p, OpenGLColor::Yellow()); Set_Point_Size(p, 5.f); }

	////boundary
	{auto p = Add_Object<OpenGLParticles<Particles<3> > >("psi_D"); Set_Visibility(p, 'd', true); Set_Color(p, OpenGLColor::Red()); }
	{auto p = Add_Object<OpenGLParticles<Particles<3> > >("psi_N"); Set_Visibility(p, 'n', true); Set_Color(p, OpenGLColor::Blue()); }
}

//////////////////////////////////////////////////////////////////////////
////solid viewers

void OpenGLViewerMesh::Toggle_Next_Frame()
{
	Base::Toggle_Next_Frame();

	////CODESAMPLE: debugging a specific vertex
	//OpenGLTriangleMesh* opengl_tri_mesh=Get_Object_By_Name_And_Type<OpenGLTriangleMesh>("triangle_mesh");
	//if(opengl_tri_mesh!=nullptr&&opengl_text_vertex!=nullptr){
	//	//Hashset<int> vtx={305,114,27,301};
	//	//opengl_text_vertex->Update_Vertex_Data(opengl_tri_mesh->mesh.Vertices(),vtx);
	//	opengl_text_vertex->Update_Vertex_Data(opengl_tri_mesh->mesh.Vertices());
	//	opengl_text_vertex->Set_Data_Refreshed();}
}

void OpenGLViewerSolid::Initialize_Data()
{
	auto p=Add_Object<OpenGLParticles<Particles<3> > >("particles");Set_Visibility(p,'x',true);Set_Color(p,OpenGLColor::Magenta());Set_Point_Size(p,4.f);
	auto psi_D=Add_Object<OpenGLParticles<Particles<3> > >("psi_D");Set_Visibility(psi_D,'x',true);Set_Color(psi_D,OpenGLColor::Gray());Set_Point_Size(psi_D,12.f);

	auto seg=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(seg,'s',true);Set_Color(seg,OpenGLColor::Red());
	auto tri=Add_Object<OpenGLTriangleMesh>("triangle_mesh");Set_Polygon_Mode(tri,PolygonMode::Fill);Set_Visibility(tri,'t',true);Set_Color(tri,OpenGLColor::Red());
	if(tri){
		tri->recomp_vtx_normal=true;
		tri->use_mat=true;
		tri->mat.mat_dif=glm::vec4(1.f,0.f,0.f,1.f);
		tri->mat.mat_spec=glm::vec4(1.f,0.f,0.f,1.f);}

	auto tet=Add_Object<OpenGLTetrahedronMesh>("tetrahedron_mesh");Set_Polygon_Mode(tet,PolygonMode::Wireframe);Set_Visibility(tet,'e',true);Set_Color(tet,OpenGLColor::Blue());

	////Set lights
	auto dir_light=OpenGLUbos::Add_Directional_Light(glm::vec3(-1.f,-.1f,-.2f));
	OpenGLUbos::Set_Ambient(glm::vec4(.1f,.1f,.1f,1.f));
	OpenGLUbos::Update_Lights_Ubo();
}

void OpenGLViewerFem::Initialize_Data()
{
	auto seg=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Polygon_Mode(seg,PolygonMode::Wireframe);
	Set_Visibility(seg,'s',true);Set_Color(seg,OpenGLColor::Red());Set_Vtx_Displacement_And_Name(seg,true,"displacement");

	auto tri=Add_Object<OpenGLTriangleMesh>("triangle_mesh");Set_Polygon_Mode(tri,PolygonMode::Wireframe);
	Set_Visibility(tri,'t',true);Set_Color(tri,OpenGLColor::Red());Set_Vtx_Displacement_And_Name(tri,true,"displacement");

	if(tri!=nullptr){
		auto strain=Add_Mesh_Tensor_Field(tri,{{"strain","",ColorType::Jet,StoreType::Cell}});Set_Visibility(strain,'4',false);
		auto stress=Add_Mesh_Tensor_Field(tri,{{"stress","",ColorType::Jet,StoreType::Cell}});Set_Visibility(stress,'5',false);
		auto von_mises=Add_Mesh_Scalar_Field(tri,{{"von_mises","",ColorType::Jet,StoreType::Cell}});Set_Visibility(von_mises,'6',true);
		auto mat=Add_Mesh_Scalar_Field(tri,{{"mat","",ColorType::Jet,StoreType::Node}});Set_Visibility(mat,'7',true);}//

	auto tet=Add_Object<OpenGLTetrahedronMesh>("tetrahedron_mesh");Set_Polygon_Mode(tet,PolygonMode::Wireframe);
	Set_Visibility(tet,'e',true);Set_Color(tet,OpenGLColor::Blue());Set_Vtx_Displacement_And_Name(tet,true,"displacement");

	if(tet!=nullptr){
		auto strain=Add_Mesh_Tensor_Field(tet,{{"strain","",ColorType::Jet,StoreType::Cell}});Set_Visibility(strain,'4',false);
		auto stress=Add_Mesh_Tensor_Field(tet,{{"stress","",ColorType::Jet,StoreType::Cell}});Set_Visibility(stress,'5',false);
		auto von_mises=Add_Mesh_Scalar_Field(tet,{{"von_mises","",ColorType::Jet,StoreType::Cell}});Set_Visibility(von_mises,'6',true);}

	auto quad=Add_Object<OpenGLQuadMesh>("quad_mesh");Set_Polygon_Mode(quad,PolygonMode::Wireframe);
	Set_Visibility(quad,'e',true);Set_Color(quad,OpenGLColor::Blue());Set_Vtx_Displacement_And_Name(quad,true,"displacement");
}

void OpenGLViewerDrone::Initialize_Data()
{
	////2D
	auto p=Add_Object<OpenGLParticles<Particles<3> > >("rotor_particles");Set_Visibility(p,'x',true);Set_Color(p,OpenGLColor::Magenta());Set_Point_Size(p,4.f);
	auto tri=Add_Object<OpenGLTriangleMesh>("body_volume_tri_mesh");Set_Visibility(tri,'t',true);Set_Color(tri,OpenGLColor::White());Set_Line_Width(tri,2.f);Set_Polygon_Mode(tri,PolygonMode::Wireframe);
	auto seg=Add_Object<OpenGLSegmentMesh>("body_surface_seg_mesh");Set_Visibility(seg,'s',true);Set_Color(seg,OpenGLColor::Blue());Set_Line_Width(seg,4.f);Set_Polygon_Mode(seg,PolygonMode::Wireframe);
	auto seg2=Add_Object<OpenGLSegmentMesh>("rotor_mesh");Set_Visibility(seg2,'s',true);Set_Color(seg2,OpenGLColor::Red());Set_Line_Width(seg2,4.f);Set_Polygon_Mode(seg2,PolygonMode::Wireframe);

	////3D
	auto tet=Add_Object<OpenGLTetrahedronMesh>("body_volume_tet_mesh");Set_Visibility(tri,'t',true);Set_Color(tri,OpenGLColor::White());Set_Line_Width(tri,2.f);Set_Polygon_Mode(tri,PolygonMode::Wireframe);
	auto tri2=Add_Object<OpenGLTriangleMesh>("body_surface_tri_mesh");Set_Visibility(seg,'s',true);Set_Color(seg,OpenGLColor::Blue());Set_Line_Width(seg,4.f);Set_Polygon_Mode(seg,PolygonMode::Wireframe);
}

//////////////////////////////////////////////////////////////////////////
////TopoOpt

void OpenGLViewerTopo::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::Gray(.4f),OpenGLColor::Gray(.2f));
}

void OpenGLViewerTopo::Initialize_Data()
{
	auto p=Add_Object<OpenGLParticles<Particles<3> > >("particles");Set_Visibility(p,'P',true);Set_Point_Size(p,3.f);

	auto grid=Add_Object<OpenGLGrid>("grid");Set_Visibility(grid,'g',true);
	//auto x=Add_Grid_Scalar_Field("grid",{{"x","",ColorType::Den,StoreType::Cell}});Set_Visibility(x,'x',true);Set_Color(x,OpenGLColor::Gray(0.f));

	if(use_2d_display){
		auto x=Add_Grid_Scalar_Field("grid",{{"x","",ColorType::Den,StoreType::Cell}});Set_Visibility(x,'x',true);Set_Color(x,OpenGLColor::Gray(0.f));}
	else{
		auto x=Add_Grid_Volume("grid",{{"x","",ColorType::Den,StoreType::Cell}});Set_Visibility(x,'x',true);}

	auto dis=Add_Grid_Vector_Field("grid",{{"displacement","1",ColorType::Jet,StoreType::Node}});Set_Visibility(dis,'v',true);
	auto von_mises=Add_Grid_Scalar_Field("grid",{{"von_mises_stress","",ColorType::Jet,StoreType::Cell}});Set_Visibility(von_mises,'6',false);
}

void OpenGLViewerCat::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::Gray(.4f),OpenGLColor::Gray(.2f));
}

void OpenGLViewerCat::Initialize_Data()
{
	////particle
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("particles");Set_Visibility(p,'f',true);
	Set_Color(p,OpenGLColor::Green());Set_Point_Size(p,4.f);
	Set_Particle_Mode(p,ParticleMode::Circle);Set_Particle_Size(p,.1f);Set_Line_Width(p,4.f);}
	
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("particles");Set_Visibility(p,'f',true);
	Set_Color(p,OpenGLColor::Blue());
	Set_Particle_Mode(p,ParticleMode::Sphere);Set_Particle_Size(p,.1f);}

	////interface
	auto c2=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(c2,'s',true);Set_Color(c2,OpenGLColor::Gray(.7f));Set_Line_Width(c2,4.f);
}

//////////////////////////////////////////////////////////////////////////
////Geometry viewers

void OpenGLViewerMesh::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::Gray(.3f),OpenGLColor::Gray(.1f));
}

void OpenGLViewerMesh::Initialize_Data()
{
	auto grid=Add_Object<OpenGLGrid>("grid");Set_Visibility(grid,'g',true);
	{auto vec=Add_Grid_Vector_Field("grid",{{"X","1",ColorType::Jet,StoreType::Node}});Set_Visibility(vec,'X',false);}
	{auto vec=Add_Grid_Vector_Field("grid",{{"force","1",ColorType::Jet,StoreType::Node}});Set_Visibility(vec,'F',false);}
	{auto tensor=Add_Grid_Tensor_Field("grid",{{"F","1",ColorType::Jet,StoreType::Node}});Set_Visibility(tensor,'T',true);}

	auto p=Add_Object<OpenGLParticles<Particles<3> > >("particles");Set_Visibility(p,'x',true);Set_Color(p,OpenGLColor::Yellow());
	Set_Particle_Velocity_Visibility(p,'v');Set_Particle_Force_Visibility(p,'f');

	auto seg=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(seg,'s',true);Set_Color(seg,OpenGLColor::Red());
	auto tri=Add_Object<OpenGLTriangleMesh>("triangle_mesh");Set_Polygon_Mode(tri,PolygonMode::Wireframe);Set_Visibility(tri,'t',true);Set_Color(tri,OpenGLColor::Red());
	auto tet=Add_Object<OpenGLTetrahedronMesh>("tetrahedron_mesh");Set_Polygon_Mode(tet,PolygonMode::Wireframe);Set_Visibility(tet,'e',true);Set_Color(tet,OpenGLColor::Blue());

	Array<std::string> texts;Array<Vector3> pos;opengl_text_vertex=Add_Interactive_Text_Array_3D(texts,pos);
}

void OpenGLViewerPoints::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::White(),OpenGLColor::White());
}

void OpenGLViewerPoints::Initialize_Data()
{
	auto p=Add_Object<OpenGLParticles<Particles<3> > >("particles");Set_Visibility(p,'x',true);
	Set_Particle_Velocity_Visibility(p,'v');Set_Particle_Force_Visibility(p,'f');
	
	Set_Color(p,OpenGLColor(0,105.f/255.f,62.f/255.f));
	Set_Particle_Mode(p,ParticleMode::Sphere);Set_Particle_Size(p,.04f);Set_Line_Width(p,4.f);//
	//Set_Polygon_Mode(p,PolygonMode::Wireframe);

	auto seg=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(seg,'s',true);
	Set_Color(seg,OpenGLColor(51.f/255.f,51.f/255.f,255.f/255.f));Set_Line_Width(seg,8.f);
}

void OpenGLViewerGeometry::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::Gray(.4f),OpenGLColor::Gray(.2f));
}

void OpenGLViewerGeometry::Initialize_Data()
{
	////grid values
	auto grid=Add_Object<OpenGLGrid>("grid");Set_Visibility(grid,'g',true);
	auto phi=Add_Grid_Scalar_Field("grid",{{"phi","",ColorType::Jet,StoreType::Cell}});Set_Visibility(phi,'l',false);

	////particle
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("particles");Set_Visibility(p,'f',true);
	Set_Color(p,OpenGLColor::Blue());Set_Point_Size(p,4.f);
	Set_Particle_Mode(p,ParticleMode::Circle);Set_Particle_Size(p,.1f);}

	////interface
	{auto c=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(c,'s',true);Set_Color(c,OpenGLColor::Red());Set_Line_Width(c,4.f);}
	{auto c=Add_Object<OpenGLTriangleMesh>("triangle_mesh");Set_Polygon_Mode(c,PolygonMode::Wireframe);Set_Color(c,OpenGLColor::Red());Set_Visibility(c,'t',true);}
	{auto c=Add_Object<OpenGLSegmentMesh>("segment_mesh_2");Set_Visibility(c,'s',true);Set_Color(c,OpenGLColor::Red());Set_Line_Width(c,4.f);}
	{auto c=Add_Object<OpenGLTriangleMesh>("triangle_mesh_2");Set_Polygon_Mode(c,PolygonMode::Wireframe);Set_Color(c,OpenGLColor::Red());Set_Visibility(c,'t',true);}
}

void OpenGLViewerCurve::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::Gray(.4f),OpenGLColor::Gray(.2f));
}

void OpenGLViewerCurve::Initialize_Data()
{
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("p1");Set_Visibility(p,'z',true);Set_Color(p,OpenGLColor::Yellow());
	auto seg=Add_Object<OpenGLSegmentMesh>("s1");Set_Visibility(seg,'x',true);Set_Color(seg,OpenGLColor::Red());}

	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("p2");Set_Visibility(p,'a',true);Set_Color(p,OpenGLColor::Green());
	auto seg=Add_Object<OpenGLSegmentMesh>("s2");Set_Visibility(seg,'s',true);Set_Color(seg,OpenGLColor::Blue());}
}

void OpenGLViewerVoronoi::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::Gray(.4f),OpenGLColor::Gray(.2f));
}

void OpenGLViewerVoronoi::Initialize_Data()
{
	auto p=Add_Object<OpenGLParticles<Particles<3> > >("particles");Set_Visibility(p,'P',true);Set_Point_Size(p,10.f);Set_Color(p,OpenGLColor::Red());
	auto p2=Add_Object<OpenGLParticles<Particles<3> > >("particles2");Set_Visibility(p2,'Q',true);Set_Point_Size(p2,10.f);Set_Color(p2,OpenGLColor::Blue());
	{auto c=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(c,'s',true);Set_Color(c,OpenGLColor::Yellow());Set_Line_Width(c,4.f);}
	{auto c=Add_Object<OpenGLTriangleMesh>("triangle_mesh");Set_Polygon_Mode(c,PolygonMode::SurfOnly);Set_Color(c,OpenGLColor::Yellow());Set_Visibility(c,'t',true);}
	auto grid=Add_Object<OpenGLGrid>("grid");Set_Visibility(grid,'g',true);

	auto vec=Add_Grid_Vector_Field("grid",{ {"grad","1",ColorType::Jet,StoreType::Cell}});Set_Visibility(vec,'v',true);Set_Scale(vec,default_scale);
	Set_Auto_Length_Norm(vec,1.);

	//if(use_2d_display){		////ATTENTION: remember to set the option "-d 2" when using the 2D mode!
	//	auto x=Add_Grid_Scalar_Field("grid",{{"rho","",ColorType::Jet,StoreType::Cell}});Set_Visibility(x,'x',true);Set_Color(x,OpenGLColor::Green());
	//	auto x2=Add_Grid_Scalar_Field("grid",{{"rho2","",ColorType::Jet,StoreType::Cell}});Set_Visibility(x2,'y' ,false);Set_Color(x2,OpenGLColor::Blue());}
	//else{
	//	auto x=Add_Grid_Volume("grid",{{"rho","",ColorType::Den,StoreType::Cell}});Set_Visibility(x,'x',true);
	//	auto x2=Add_Grid_Volume("grid",{{"rho2","",ColorType::Den,StoreType::Cell}});Set_Visibility(x2,'y',true);}

	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("psi_D");Set_Visibility(p,'D',true);Set_Color(p,OpenGLColor::Red());Set_Scale(p,default_scale);}
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("psi_N");Set_Visibility(p,'N',true);Set_Color(p,OpenGLColor::Blue());Set_Particle_Force_Visibility(p,'F',true);Set_Scale(p,default_scale);}

	auto x=Add_Grid_Scalar_Field("grid",{{"rho","",ColorType::Jet,StoreType::Cell}});Set_Visibility(x,'x',true);Set_Color(x,OpenGLColor::Green());
	auto x2=Add_Grid_Scalar_Field("grid",{{"rho2","",ColorType::Jet,StoreType::Cell}});Set_Visibility(x2,'y',false);Set_Color(x2,OpenGLColor::Blue());

	auto dis=Add_Grid_Vector_Field("grid",{{"displacement","1",ColorType::Jet,StoreType::Node}});Set_Visibility(dis,'d',true);Set_Scale(dis,default_scale);
	auto phi=Add_Grid_Scalar_Field("grid",{{"phi","",ColorType::Jet,StoreType::Cell}});Set_Visibility(phi,'L',false);
	auto active=Add_Grid_Scalar_Field("grid",{{"active","",ColorType::Jet,StoreType::Cell}});Set_Visibility(active,'a',false);
}

void OpenGLViewerVortex::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::Gray(.4f),OpenGLColor::Gray(.2f));
}

void OpenGLViewerVortex::Initialize_Data()
{
	////grid values
	auto grid=Add_Object<OpenGLGrid>("grid");Set_Visibility(grid,'g',true);
	{auto vec=Add_Grid_Vector_Field("grid",{{"velocity","1",ColorType::Jet,StoreType::Cell}});Set_Visibility(vec,'v',true);Set_Scale(vec,default_scale);}
	{auto vec=Add_Grid_Vector_Field("grid",{ {"impulse","",ColorType::Jet,StoreType::Cell}});Set_Visibility(vec,'a',true);Set_Scale(vec,default_scale);}
	{auto psi=Add_Grid_Vector_Field("grid",{{"psi","",ColorType::Jet,StoreType::Cell}});Set_Visibility(psi,'S',true);Set_Scale(psi,default_scale);}

	//auto vor=Add_Grid_Scalar_Field("grid",{{"vor","",ColorType::Jet,StoreType::Node}});Set_Visibility(vor,'o',false);
	//auto vor3d=Add_Grid_Vector_Field("grid",{{"vor3d","",ColorType::Jet,StoreType::Node}});Set_Visibility(vor3d,'x',false);
	
	//auto s=Add_Grid_Vector_Field("grid",{{"s","",ColorType::Jet,StoreType::Cell}});Set_Visibility(s,'S',true);
	//auto pressure=Add_Grid_Scalar_Field("grid",{{"pressure","",ColorType::Jet,StoreType::Cell}});Set_Visibility(pressure,'i',false);

	////particle
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("particles");Set_Visibility(p,'f',true);
	Set_Color(p,OpenGLColor::Blue());Set_Point_Size(p,2.f);Set_Particle_Mode(p,ParticleMode::Dot);
	Set_Particle_Velocity_Visibility(p,'k');
	Set_Particle_Force_Visibility(p,'n');}

	////boundary
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("psi_D");Set_Visibility(p,'d',true);Set_Color(p,OpenGLColor::Red());}
	{auto p=Add_Object<OpenGLParticles<Particles<3> > >("psi_N");Set_Visibility(p,'N',true);Set_Color(p,OpenGLColor::Blue());}
}

void OpenGLViewerImpulse::Initialize()
{
	Base::Initialize();
	//Set_Background_Color(OpenGLColor::White(),OpenGLColor::White());
	Set_Background_Color(OpenGLColor::Gray(.4f),OpenGLColor::Gray(.2f));
}

void OpenGLViewerImpulse::Initialize_Data()
{
	switch(test){
	case 1:{	////standard view
		////grid and fields
		auto grid=Add_Object<OpenGLGrid>("grid");Set_Visibility(grid,'g',true);
		if(use_2d_display){auto density=Add_Grid_Scalar_Field("grid",{{"density","",ColorType::Den,StoreType::Cell}});Set_Visibility(density,'d',true);}
		else{auto density=Add_Grid_Volume("grid",{{"density","",ColorType::Den,StoreType::Cell}});Set_Visibility(density,'d',true);}
		if(!minimal_mode){
			////fields
			{auto vec=Add_Grid_Vector_Field("grid",{ {"velocity","",ColorType::Jet,StoreType::Cell}});Set_Visibility(vec,'v',true);Set_Scale(vec,default_scale);}
			{auto vec=Add_Grid_Vector_Field("grid",{ {"impulse","",ColorType::Jet,StoreType::Cell}});Set_Visibility(vec,'a',true);Set_Scale(vec,default_scale);}	
			{auto vor=Add_Grid_Scalar_Field("grid",{{"vorticity","",ColorType::Jet,StoreType::Cell}});Set_Visibility(vor,'o',false);}
			{auto solid_phi=Add_Grid_Scalar_Field("grid",{{"solid_phi","",ColorType::Jet,StoreType::Cell}});Set_Visibility(solid_phi,'S',false);}
			{auto fluid_phi=Add_Grid_Scalar_Field("grid",{{"fluid_phi","",ColorType::Jet,StoreType::Cell}});Set_Visibility(fluid_phi,'F',false);}
			{auto p=Add_Grid_Scalar_Field("grid",{{"p","",ColorType::Jet,StoreType::Cell}});Set_Visibility(p,'P',false);}
			{auto q=Add_Grid_Scalar_Field("grid",{{"q","",ColorType::Jet,StoreType::Cell}});Set_Visibility(q,'Q',false);}
			{auto vec=Add_Grid_Vector_Field("grid",{ {"grad_q","",ColorType::Jet,StoreType::Cell}});Set_Visibility(vec,'X',false);Set_Scale(vec,default_scale);}
			{auto vec=Add_Grid_Vector_Field("grid",{ {"grad_p","",ColorType::Jet,StoreType::Cell}});Set_Visibility(vec,'Z',false);Set_Scale(vec,default_scale);}	
			{auto type=Add_Grid_Scalar_Field("grid",{{"fluid_type","",ColorType::Jet,StoreType::Cell}});Set_Visibility(type,'T',false);}
			{auto psi=Add_Grid_Vector_Field("grid",{{"psi","",ColorType::Jet,StoreType::Cell}});Set_Visibility(psi,'S',true);Set_Scale(psi,default_scale);}

			{auto p=Add_Object<OpenGLParticles<Particles<3> > >("solid_particles");Set_Visibility(p,'f',true);
			Set_Color(p,OpenGLColor::Black());Set_Point_Size(p,12.f);/*Set_Particle_Velocity_Visibility(p,'v');Set_Particle_Force_Visibility(p,'F');*/}

			////mesh for levelset
			auto c2=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(c2,'s',true);Set_Color(c2,OpenGLColor::Red());Set_Line_Width(c2,4.f);
			auto c3=Add_Object<OpenGLTriangleMesh>("triangle_mesh");Set_Polygon_Mode(c3,PolygonMode::Wireframe);Set_Shading_Mode(c3,ShadingMode::None);
			Set_Color(c3,OpenGLColor::Green(.5f));Set_Visibility(c3,'t',true);Set_Line_Width(c3,4.f);
	
			////boundary conditions
			{auto p=Add_Object<OpenGLParticles<Particles<3> > >("psi_D");Set_Visibility(p,'D',true);Set_Color(p,OpenGLColor::Red());}
			{auto p=Add_Object<OpenGLParticles<Particles<3> > >("psi_N");Set_Visibility(p,'N',true);Set_Color(p,OpenGLColor::Blue());}

			////particle
			{auto p=Add_Object<OpenGLParticles<Particles<3> > >("particles");Set_Visibility(p,'f',true);
			Set_Color(p,OpenGLColor::Red());Set_Point_Size(p,8.f);Set_Particle_Velocity_Visibility(p,'v');Set_Particle_Force_Visibility(p,'F');}

			{auto p=Add_Object<OpenGLTrackerPoints>("tracker_points");Set_Visibility(p,'R',true);Set_Color(p,OpenGLColor::Yellow());Set_Point_Size(p,4.f);}
			{auto p=Add_Object<OpenGLTrackerPoints>("tracker_points_2");Set_Visibility(p,'U',true);Set_Color(p,OpenGLColor::Red());Set_Point_Size(p,6.f);}
			//{auto p=Add_Object<OpenGLTrackerVectors>("point_force");Set_Visibility(p,'F',true);Set_Color(p,OpenGLColor::Green());Set_Scale(p,default_scale);}
			{auto p=Add_Object<OpenGLTrackerVectors>("point_velocity");Set_Visibility(p,'V',true);Set_Color(p,OpenGLColor::Blue());Set_Scale(p,default_scale);}
			{auto p=Add_Object<OpenGLTrackerVectors>("point_impulse");Set_Visibility(p,'I',true);Set_Color(p,OpenGLColor::Red());Set_Scale(p,default_scale);}

			////faces
			{auto p=Add_Object<OpenGLTrackerVectors>("velocity_face");Set_Visibility(p,'V',false);Set_Color(p,OpenGLColor::Blue());Set_Scale(p,default_scale);}
			{auto p=Add_Object<OpenGLTrackerVectors>("impulse_face");Set_Visibility(p,'A',false);Set_Color(p,OpenGLColor::Red());Set_Scale(p,default_scale);}}	
	}break;
	case 2:{	////particle view only
		auto p=Add_Object<OpenGLTrackerPoints>("marker_particles");Set_Visibility(p,'x',true);Set_Color(p,OpenGLColor::Yellow());Set_Point_Size(p,1.f);
	}break;
	}
}

void OpenGLViewerImpulseDemo::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::White(),OpenGLColor::White());
}

void OpenGLViewerImpulseDemo::Initialize_Data()
{
	if(use_2d_display){
		auto c2=Add_Object<OpenGLSegmentMesh>("segment_mesh");Set_Visibility(c2,'s',true);Set_Color(c2,OpenGLColor::Black());Set_Line_Width(c2,6.f);
		auto c21 = Add_Object<OpenGLSegmentMesh>("woodfish_segment_mesh"); Set_Visibility(c21, 'S', true); Set_Color(c21, OpenGLColor::Black()); Set_Line_Width(c21, 4.f);
		auto c22 = Add_Object<OpenGLTriangleMesh>("woodfish_surface_mesh"); Set_Polygon_Mode(c22, PolygonMode::Fill); Set_Shading_Mode(c22, ShadingMode::None);
		Set_Color(c22, OpenGLColor::Yellow(.8f)); Set_Visibility(c22, 'T', true); Set_Line_Width(c22, 4.f);
	}else{
		auto c3=Add_Object<OpenGLTriangleMesh>("triangle_mesh");Set_Polygon_Mode(c3,PolygonMode::Wireframe);Set_Shading_Mode(c3,ShadingMode::None);
		Set_Color(c3,OpenGLColor::Green(.5f));Set_Visibility(c3,'t',true);Set_Line_Width(c3,4.f);
		auto c31 = Add_Object<OpenGLTriangleMesh>("woodfish_triangle_mesh"); Set_Polygon_Mode(c31, PolygonMode::Wireframe); Set_Shading_Mode(c31, ShadingMode::None);
		Set_Color(c31, OpenGLColor::Green(.5f)); Set_Visibility(c31, 'T', true); Set_Line_Width(c31, 4.f);
	}

	//{auto p=Add_Object<OpenGLTrackerPoints>("tracker_points");Set_Visibility(p,'R',true);Set_Color(p,OpenGLColor::Black());Set_Point_Size(p,2.f);}
	{auto p=Add_Object<OpenGLTrackerPoints>("colored_tracker_points");Set_Visibility(p,'R',true);Set_Color(p,OpenGLColor::Black());Set_Point_Size(p,4.f);
	if(p)p->Use_Point_Color();}

	{auto p=Add_Object<OpenGLTrackerVectors>("point_velocity");Set_Visibility(p,'V',true);Set_Color(p,OpenGLColor::Green(.5f));Set_Scale(p,default_scale);}
	{auto p=Add_Object<OpenGLTrackerVectors>("point_impulse");Set_Visibility(p,'I',true);Set_Color(p,OpenGLColor::Blue());Set_Scale(p,default_scale);}

	auto opengl_segments=Add_Interactive_Object<OpenGLSegments>();
	opengl_segments->vertices.push_back(Vector3(-1.,0.001,0.));
	opengl_segments->vertices.push_back(Vector3(3.,0.001,0.));
	opengl_segments->color=OpenGLColor(0.3f,0.3f,0.3f,1.f);
	opengl_segments->line_width=8.f;
	opengl_segments->Set_Data_Refreshed();
	opengl_segments->Initialize();
}

void OpenGLViewerElasticityDEC::Initialize()
{
	Base::Initialize();
	Set_Background_Color(OpenGLColor::Gray(.4f), OpenGLColor::Gray(.2f));
}

void OpenGLViewerElasticityDEC::Initialize_Data()
{
	////grid values
	auto grid = Add_Object<OpenGLGrid>("grid"); Set_Visibility(grid, 'g', true);

	auto dis = Add_Grid_Vector_Field("grid", { {"displacement","1",ColorType::Jet,StoreType::Node} }); Set_Visibility(dis, 'D', true);
	//Set_Auto_Length_Norm(dis, 1.);

	auto vec = Add_Grid_Vector_Field("grid", { {"velocity","1",ColorType::Jet,StoreType::Node} }); Set_Visibility(vec, 'v', false);

	auto f = Add_Grid_Vector_Field("grid", { {"force","1",ColorType::Jet,StoreType::Node} }); Set_Visibility(f, 'F', false);

	auto x = Add_Grid_Scalar_Field("grid", { {"x","",ColorType::Den,StoreType::Cell} }); Set_Visibility(x, 'x', false); Set_Color(x, OpenGLColor::Red());

	auto y = Add_Grid_Volume("grid", { {"x","",ColorType::Den,StoreType::Cell} }); Set_Visibility(y, 'y', false);
}
