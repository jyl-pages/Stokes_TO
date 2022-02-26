#ifndef __SketchPadDriver_h__
#define __SketchPadDriver_h__
#include <random>
#include "Common.h"
#include "Driver.h"
#include "OpenGLMarkerObjects.h"
#include "OpenGLMesh.h"
#include "OpenGLCommon.h"
#include "OpenGLWindow.h"
#include "OpenGLViewer.h"
#include "OpenGLParticles.h"
#include "MeshAdvFunc.h"
#include "Triangulation2D.h"

template<int d> class SketchPadDriver : public Driver, public OpenGLViewer
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;using Base=Driver;
	OpenGLTriangleMesh* opengl_triangle_mesh=nullptr;
	OpenGLParticles<Particles<3> >* opengl_particles=nullptr;
	TriangleMesh<2> tri_mesh;
	
public:
	virtual void Initialize()
	{
		////viewer initialization, initialize visualization data
		OpenGLViewer::Initialize();
	}

	////synchronize simulation data to visualization data
	virtual void Initialize_Data()
	{
		//MeshFunc::Initialize_Circle_Mesh((real)1,&tri_mesh,32);
		MeshFunc::Initialize_Rectangle_Mesh((real)2,(real)1,Vector2::Zero(),&tri_mesh,(real).2);

		opengl_triangle_mesh=Add_Interactive_Object<OpenGLTriangleMesh>();
		Set_Polygon_Mode(opengl_triangle_mesh,PolygonMode::Wireframe);

		opengl_particles=Add_Interactive_Object<OpenGLParticles<Particles<3> > >();

		Sync_Simulation_And_Visualization_Data();
	}

	void Sync_Simulation_And_Visualization_Data()
	{
		auto& mesh=opengl_triangle_mesh->mesh;
		mesh.Vertices().clear();
		mesh.Elements().clear();

		for(int i=0;i<tri_mesh.Vertices().size();i++){
			mesh.Vertices().push_back(V3(tri_mesh.Vertices()[i]));}
		mesh.Elements()=tri_mesh.Elements();

		opengl_triangle_mesh->Set_Data_Refreshed();
		opengl_triangle_mesh->Initialize();

		auto& particles=opengl_particles->particles;
		particles.Resize((int)tri_mesh.Vertices().size());
		for(int i=0;i<tri_mesh.Vertices().size();i++){
			particles.X(i)=V3(tri_mesh.Vertices()[i]);}

		opengl_particles->Set_Data_Refreshed();
		opengl_particles->Initialize();
	}

	////update simulation and visualization for each time step
	virtual void Toggle_Next_Frame()
	{
		Sync_Simulation_And_Visualization_Data();
		OpenGLViewer::Toggle_Next_Frame();
	}

	virtual void Run()
	{
		OpenGLViewer::Run();
	}

	////User interaction
	virtual bool Mouse_Click(int left,int right,int mid,int x,int y,int w,int h)
	{
		if(left!=1){return false;}
		Vector3f win_pos=opengl_window->Project(Vector3f::Zero());
		Vector3f pos=opengl_window->Unproject(Vector3f((float)x,(float)y,win_pos[2]));
		VectorD p_pos;for(int i=0;i<d;i++)p_pos[i]=(real)pos[i];
		real r=.1*static_cast<float>(rand()%1000)/1000.+.15;
		std::cout<<"click pos: "<<p_pos.transpose()<<std::endl;
		real r_sq=pow((real).2,2);
		for(int i=0;i<tri_mesh.Vertices().size();i++){
			if((p_pos-tri_mesh.Vertices()[i]).norm()<r_sq)std::cout<<"select "<<i<<std::endl;}
		return true;
	}


protected:
	////Helper function to convert a vector to 3d, for c++ template
	Vector3 V3(const Vector2& v2){return Vector3(v2[0],v2[1],.0);}
	Vector3 V3(const Vector3& v3){return v3;}
};
#endif
