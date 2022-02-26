//////////////////////////////////////////////////////////////////////////
// Opengl viewer
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLViewer_h__
#define __OpenGLViewer_h__
#include <iostream>
#include <memory>
#include "Hashtable.h"
#include "Common.h"
#include "OpenGLObject.h"
#include "OpenGLAabb.h"

////Forward declaration
class OpenGLWindow;
class OpenGLGrid;
template<class T,int data_size=4> class OpenGLGridField;
template<class T> class OpenGLGridVectors;
template<class T> class OpenGLGridTensors;
template<class T> class OpenGLGridHeightField;
template<int d> class TriangleMesh;
template<int d> class TetrahedronMesh;
class OpenGLSegmentMesh;
class OpenGLTriangleMesh;
class OpenGLTetrahedronMesh;
template<class T,class T_MESH,int data_size=4> class OpenGLMeshField;
template<class T,class T_MESH> class OpenGLMeshVectors;
template<class T,class T_MESH> class OpenGLMeshTensors;
class OpenGLUICommand;
class OpenGLUICurvePlot;
class OpenGLBackground;
class OpenGLText;
class OpenGLText3D;
class OpenGLTextArray3D;
template<class T_PARTICLES> class OpenGLParticles;
template<int d,typename T> class Particles;
class OpenGLSphere;
class OpenGLArrow;

class OpenGLViewer
{
public:
	std::string output_dir;
	std::string config_file_name;
	std::shared_ptr<OpenGLWindow> opengl_window=nullptr;
	int first_frame=0,last_frame=-1,frame=0;
	int every=1;
	bool draw_bk=true;
	bool draw_axes=false;
	bool draw_box=false;
	bool use_ui=true;
	bool play=false;
	bool use_2d_display=false;
	bool verbose=false;
	bool minimal_mode=false;
	HashtableMultiValue<uchar,std::string> key_data_hashtable;
	std::shared_ptr<OpenGLUICommand> ui_command=nullptr;
	real default_scale=(real)1;
	int test=1;		////specifying different test modes

	//////////////////////////////////////////////////////////////////////////
	////Initialization and run
	virtual void Initialize();
	virtual void Initialize_Data(){}
	virtual void Run();
	virtual void Initialize_Common_Data();
	virtual void Initialize_Camera();
	////Customized mouse func, the same as OpenGLObject
	virtual bool Mouse_Drag(int x,int y,int w,int h){return false;}
	virtual bool Mouse_Click(int left,int right,int mid,int x,int y,int w,int h){return false;}

	//////////////////////////////////////////////////////////////////////////
	////Animation
	void Update_Frame();

	//////////////////////////////////////////////////////////////////////////
	////UI
	virtual void Initialize_UI();
	virtual void Finish();
	virtual void Toggle_Command(const std::string cmd);		////supporting multiple commands separated by ';'
	void Toggle_Single_Command(const std::string cmd);		////the previous Toggle_Command for a single command

	//////////////////////////////////////////////////////////////////////////
	////Add objects
	template<class T_OBJECT> T_OBJECT* Add_Object(const std::string object_name,const Array<OpenGLData> data=Array<OpenGLData>())
	{
		T_OBJECT* opengl_object=nullptr;
		if(Initialize_From_File(output_dir,opengl_object,object_name,data,first_frame)){
			Add_OpenGL_Object(opengl_object);
			if(verbose){std::cout<<"Add opengl object: "<<object_name;}}
		return opengl_object;
	}

	template<class T_OBJECT> T_OBJECT* Add_Object(const char* object_name,const Array<OpenGLData> data=Array<OpenGLData>())
	{return Add_Object<T_OBJECT>(std::string(object_name),data);}

	template<class T_OBJECT> T_OBJECT* Add_Object(T_OBJECT* opengl_object,bool init=true,bool interactive=false)
	{
		if(init)opengl_object->Initialize();
		opengl_object->interactive=interactive;
		Add_OpenGL_Object(opengl_object);
		return opengl_object;
	}

	template<class T_OBJECT> T_OBJECT* Add_Object(bool init=true,bool interactive=false)
	{
		T_OBJECT* opengl_object=new T_OBJECT();
		return Add_Object<T_OBJECT>(opengl_object,init,interactive);
	}

	template<class T_OBJECT> T_OBJECT* Add_Interactive_Object(bool init=false)
	{
		return Add_Object<T_OBJECT>(init,true);
	}

	void Add_OpenGL_Object(OpenGLObject* object);

	OpenGLGridField<real>* Add_Grid_Scalar_Field(const std::string& grid_name,const Array<OpenGLData> data=Array<OpenGLData>());
	OpenGLGridVectors<real>* Add_Grid_Vector_Field(const std::string& grid_name,const Array<OpenGLData> data=Array<OpenGLData>());
	OpenGLGridTensors<real>* Add_Grid_Tensor_Field(const std::string& grid_name,const Array<OpenGLData> data=Array<OpenGLData>());
	OpenGLGridHeightField<real>* Add_Grid_Height_Field(const std::string& grid_name,const Array<OpenGLData> data=Array<OpenGLData>());

	OpenGLMeshField<real,TriangleMesh<3> >* Add_Mesh_Scalar_Field(OpenGLTriangleMesh* opengl_mesh,const Array<OpenGLData> data=Array<OpenGLData>());
	OpenGLMeshField<real,TetrahedronMesh<3> >* Add_Mesh_Scalar_Field(OpenGLTetrahedronMesh* opengl_mesh,const Array<OpenGLData> data=Array<OpenGLData>());
	OpenGLMeshVectors<real,TriangleMesh<3> >* Add_Mesh_Vector_Field(OpenGLTriangleMesh* opengl_mesh,const Array<OpenGLData> data=Array<OpenGLData>());
	OpenGLMeshVectors<real,TetrahedronMesh<3> >* Add_Mesh_Vector_Field(OpenGLTetrahedronMesh* opengl_mesh,const Array<OpenGLData> data=Array<OpenGLData>());
	OpenGLMeshTensors<real,TriangleMesh<3> >* Add_Mesh_Tensor_Field(OpenGLTriangleMesh* opengl_mesh,const Array<OpenGLData> data=Array<OpenGLData>());
	OpenGLMeshTensors<real,TetrahedronMesh<3> >* Add_Mesh_Tensor_Field(OpenGLTetrahedronMesh* opengl_mesh,const Array<OpenGLData> data=Array<OpenGLData>());

	OpenGLVolume<real>* Add_Grid_Volume(const std::string& grid_name,const Array<OpenGLData> data=Array<OpenGLData>());

	template<class T_OBJECT> bool Initialize_From_File(const std::string& output_dir,T_OBJECT* & opengl_object,
		const std::string& object_name,const Array<OpenGLData>& data,const int frame=0)
	{
		if(OpenGLObject::Object_File_Exists(output_dir,frame,object_name)){
			opengl_object=new T_OBJECT();
			opengl_object->output_dir=output_dir;
			opengl_object->name=object_name;
			opengl_object->data=data;
			opengl_object->Refresh(frame);return true;}
		return false;
	}

	//////////////////////////////////////////////////////////////////////////
	////Interactive object creation and access
	OpenGLText* Add_Interactive_Text_2D(const std::string& text,const Vector2& pos);
	OpenGLText3D* Add_Interactive_Text_3D(const std::string& text,const Vector3& pos);
	OpenGLTextArray3D* Add_Interactive_Text_Array_3D(const Array<std::string>& text,const Array<Vector3>& pos);
	OpenGLTriangleMesh* Add_Interactive_Triangle_Mesh();
	OpenGLParticles<Particles<3,real> >* Add_Interactive_Particles();
	OpenGLSphere* Add_Interactive_Sphere(const Vector3& center=Vector3::Zero(),const real r=(real)1);
	OpenGLArrow* Add_Interactive_Arrow(const Vector3& start=Vector3::Zero(),const Vector3& end=Vector3::Unit(1));

	OpenGLObject* Get_Object_By_Name(const std::string& name);
	template<class T_OBJECT> T_OBJECT* Get_Object_By_Name_And_Type(const std::string& name)
	{OpenGLObject* obj=Get_Object_By_Name(name);if(obj==nullptr)return nullptr;T_OBJECT* obj_with_type=dynamic_cast<T_OBJECT*>(obj);return obj_with_type;}

	OpenGLBackground* Get_Object_Background();
	OpenGLTriangleMesh* Get_Object_Triangle_Mesh_By_Name(const std::string& name);

	void Set_Background_Color(const OpenGLColor& c1,const OpenGLColor& c2);

	//////////////////////////////////////////////////////////////////////////
	////Set viewer and object properties
	template<class T_OBJECT> void Set_Visibility(T_OBJECT* obj,const uchar key,const bool visible=true,const std::string draw="")
	{if(obj==nullptr)return;obj->visible=visible;Bind_Draw_Callback_Key(key,obj,draw);}
	
	template<class T_OBJECT> void Set_Color(T_OBJECT* obj,const OpenGLColor& color)
	{if(obj==nullptr)return;obj->Set_Color(color);}

	template<class T_OBJECT> void Set_Multi_Color(T_OBJECT* obj,const Array<OpenGLColor>& multi_color)
	{if(obj==nullptr)return;obj->multi_color=multi_color;}

	template<class T_OBJECT> void Set_Multi_Color(T_OBJECT* obj,const Array<float>& multi_color)
	{if(obj==nullptr)return;int n=(int)multi_color.size()/4;obj->multi_color.clear();
	for(int i=0;i<n;i++)obj->multi_color.push_back(OpenGLColor(multi_color[i*4],multi_color[i*4+1],multi_color[i*4+2],multi_color[i*4+3]));}

	template<class T_OBJECT> void Set_Alpha(T_OBJECT* obj,const GLfloat alpha)
	{if(obj==nullptr)return;obj->alpha=alpha;}

	template<class T_OBJECT> void Set_Data_Color(T_OBJECT* obj,const int data_idx,const OpenGLColor& color)
	{if(obj==nullptr||data_idx>=obj->data.size())return;for(int i=0;i<4;i++)obj->data[data_idx].color[i]=color.rgba[i];}

	template<class T_OBJECT> void Set_Scale(T_OBJECT* obj,const real scale)
	{if(obj==nullptr)return;obj->scale=scale;}

	template<class T_OBJECT> void Set_Polygon_Mode(T_OBJECT* obj,PolygonMode polygon_mode)
	{if(obj==nullptr)return;obj->polygon_mode=polygon_mode;}

	template<class T_OBJECT> void Set_Shading_Mode(T_OBJECT* obj,ShadingMode shading_mode)
	{if(obj==nullptr)return;obj->Set_Shading_Mode(shading_mode);}

	template<class T_OBJECT> void Set_Particle_Mode(T_OBJECT* obj,ParticleMode particle_mode)
	{if(obj==nullptr)return;obj->Set_Particle_Mode(particle_mode);}
	
	template<class T_OBJECT> void Set_Line_Width(T_OBJECT* obj,const GLfloat line_width)
	{if(obj==nullptr)return;obj->Set_Line_Width(line_width);}

	template<class T_OBJECT> void Set_Point_Size(T_OBJECT* obj,const GLfloat point_size)
	{if(obj==nullptr)return;obj->Set_Point_Size(point_size);}

	template<class T_OBJECT> void Set_Particle_Size(T_OBJECT* obj,const GLfloat particle_size)
	{if(obj==nullptr)return;obj->Set_Particle_Size((real)particle_size);}

	template<class T_OBJECT> void Set_Subgrid_Slice(T_OBJECT* obj,const int axis,const int start,const int end=-1)
	{if(obj==nullptr)return;obj->Set_Subgrid_Slice(axis,start,end);}

	template<class T_OBJECT> void Set_Auto_Length_Norm(T_OBJECT* obj,const real factor=(real)1)
	{if(obj==nullptr)return;obj->line_norm=obj->mac_grid.grid.dx*factor;obj->Normalize_Data();}

	void Set_Offscreen_Output_Dir(const std::string _offscreen_output);

	template<class T_OBJECT> void Set_Use_Vtx_Displacement(T_OBJECT* obj,const bool use_vtx_displacement=true,const std::string displacement_name="displacement")
	{if(obj==nullptr)return;obj->use_vtx_displacement=use_vtx_displacement;obj->displacement_name=displacement_name;}

	template<class T_OBJECT> void Set_Displacement_Name(T_OBJECT* obj,const std::string displacement_name="displacement")
	{if(obj==nullptr)return;obj->displacement_name=displacement_name;}

	template<class T_OBJECT> void Set_Vtx_Displacement_And_Name(T_OBJECT* obj,const bool use_vtx_displacement,const std::string displacement_name)
	{if(obj==nullptr)return;obj->use_vtx_displacement=use_vtx_displacement;obj->displacement_name=displacement_name;}

	template<class T_OBJECT> void Set_Particle_Velocity_Visibility(T_OBJECT* obj,const uchar key,const bool visible=true,const std::string draw="draw particle v")
	{if(obj==nullptr)return;obj->Set_Visibility_For_Velocity_Field();Bind_Callback_Key(key,&obj->Toggle_Draw_Velocity_Func,draw);}

	template<class T_OBJECT> void Set_Particle_Force_Visibility(T_OBJECT* obj,const uchar key,const bool visible=true,const std::string draw="draw particle f")
	{if(obj==nullptr)return;obj->Set_Visibility_For_Force_Field();Bind_Callback_Key(key,&obj->Toggle_Draw_Force_Func,draw);}

	//////////////////////////////////////////////////////////////////////////
	////Basic callbacks
	void Print_Keyboard_Callbacks();
	// Generates a std::function called Print_Keyboard_Callbacks_Func
	// Similar for others
	Define_Function_Object(OpenGLViewer,Print_Keyboard_Callbacks); 

	void Object_Info_Debug();
	Define_Function_Object(OpenGLViewer, Object_Info_Debug);

	virtual void Toggle_Increase_Scale();
	Define_Function_Object(OpenGLViewer,Toggle_Increase_Scale);

	virtual void Toggle_Decrease_Scale();
	Define_Function_Object(OpenGLViewer,Toggle_Decrease_Scale);

	virtual void Toggle_Increase_Feature();
	Define_Function_Object(OpenGLViewer, Toggle_Increase_Feature);

	virtual void Toggle_Decrease_Feature();
	Define_Function_Object(OpenGLViewer, Toggle_Decrease_Feature);

	virtual void Toggle_Normalize_Data();
	Define_Function_Object(OpenGLViewer,Toggle_Normalize_Data);

	virtual void Toggle_Draw_Dis();
	Define_Function_Object(OpenGLViewer,Toggle_Draw_Dis);

	virtual void Toggle_Next_Frame();
	Define_Function_Object(OpenGLViewer,Toggle_Next_Frame);

	virtual void Toggle_Prev_Frame();
	Define_Function_Object(OpenGLViewer,Toggle_Prev_Frame);

	virtual void Toggle_First_Frame();
	Define_Function_Object(OpenGLViewer,Toggle_First_Frame);

	virtual void Toggle_Play();
	Define_Function_Object(OpenGLViewer,Toggle_Play);

	virtual void Initialize_Common_Callback_Keys();

	////Customized callback function: Toggle_Func_Idx_{0-3}
	Define_Function_Object_Toggle_Func_Idx(OpenGLViewer,0);
	Define_Function_Object_Toggle_Func_Idx(OpenGLViewer,1);
	Define_Function_Object_Toggle_Func_Idx(OpenGLViewer,2);
	Define_Function_Object_Toggle_Func_Idx(OpenGLViewer,3);

	void Bind_Callback_Key(const uchar key,std::function<void(void)>* callback,const std::string& discription);
	
	void Bind_Callback_Keys(const Array<OpenGLData>& data,Array<std::function<void(void)>*> data_idx_callbacks,int start_idx=0);

	template<class T_OBJECT> void Bind_Draw_Callback_Key(const uchar key,T_OBJECT* obj,const std::string _draw="")
	{if(obj==nullptr)return;std::string draw=_draw;if(draw=="")draw="draw "+obj->name;Bind_Callback_Key(key,&obj->Toggle_Draw_Func,draw);}

	template<class T_OBJECT> void Bind_Func_Idx_0_Callback_Key(const uchar key,T_OBJECT* obj,const std::string _func="")
	{if(obj==nullptr)return;std::string func=_func;if(func=="")func="func "+obj->name;Bind_Callback_Key(key,&obj->Toggle_Func_Idx_0_Func,func);}

	//template<class T_OBJECT> OpenGLColorBar* Bind_Color_Bar(T_OBJECT* obj,const uchar key,const Vector2 bar_pos=Vector2((real)100,(real)100))
	//{
	//	if(obj==nullptr)return nullptr;
	//	OpenGLColorBar* b=new OpenGLColorBar(bar_pos);
	//	b->Set_Data_Pointer(&(obj->color_mapper),&(obj->data_idx));Add_Object(b);
	//	obj->binded_objects.push_back(b);
	//	Set_Visibility(b,key,obj->visible,"draw color bar");
	//	return b;
	//}

protected:
};
#endif