//////////////////////////////////////////////////////////////////////////
// Opengl viewer
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <memory>
#include <sstream>
#include "Mesh.h"
#include "File.h"
#include "Particles.h"
#include "OpenGLWindow.h"
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
#include "OpenGLScreenObjects.h"
#include "OpenGLMarkerObjects.h"
#include "OpenGLUi.h"
#include "OpenGLParticles.h"

//////////////////////////////////////////////////////////////////////////
////Initialization and run

void OpenGLViewer::Initialize()
{
	opengl_window=std::make_shared<OpenGLWindow>(); // Create OpenGLWindow
    opengl_window->Init(); // Initialize the window object
	opengl_window->opengl_viewer.reset(this); // Set the opengl_viewer pointer in it
	Initialize_Common_Callback_Keys(); // Like: 'p' for play and 'q' for quit
	Initialize_Common_Data(); // Add background, bounding box, axes, therefore product frame 0
	Initialize_Data(); // Depending on specific viewer. Like: velocity field , meshes...
	Initialize_Camera(); // Nothing happens here
	Initialize_UI(); // Just initialize UI
	Update_Frame(); // Read a frame from folder and send them to output
	//opengl_window->Update_Data_To_Render();
	Print_Keyboard_Callbacks();
}

void OpenGLViewer::Run()
{
	opengl_window->Run();
}

void OpenGLViewer::Add_OpenGL_Object(OpenGLObject* object)
{
	opengl_window->Add_Object(object);
}

void OpenGLViewer::Initialize_Common_Data()
{
	opengl_window->use_2d_display=use_2d_display;
	if(draw_bk)Add_Interactive_Object<OpenGLBackground>(/*init*/true);
	if(draw_axes){
		auto axes=Add_Interactive_Object<OpenGLAxes>();
		axes->use_2d_display=use_2d_display;axes->Initialize();}
	if(draw_box){
		auto box=Add_Interactive_Object<OpenGLSquare>();
		box->Initialize();}

	std::string file_name=output_dir+"/0/last_frame.txt";
	if(File::File_Exists(file_name))
		File::Read_Text_From_File(file_name,last_frame);
	if(verbose)std::cout<<"Read last frame: "<<last_frame<<std::endl;
}

void OpenGLViewer::Initialize_Camera()
{
	//opengl_window->splined_camera.Initialize();
}

//////////////////////////////////////////////////////////////////////////
////Animation

void OpenGLViewer::Update_Frame()
{
	//Read non-interactive objects from frame folder
	//Then update all data to render
	for (auto& obj : opengl_window->object_list) {
		if (!obj->interactive)obj->Refresh(frame);
		obj->Update_Data_To_Render();
	}

	//Update frame number
	opengl_window->texts["frame"]="Frame: "+std::to_string(frame);

	if(opengl_window->display_offscreen){
		opengl_window->frame_offscreen=frame;
		opengl_window->frame_offscreen_rendered=-1;}

	opengl_window->Update_Color_Code_Offset();
}

//////////////////////////////////////////////////////////////////////////
////UI

void OpenGLViewer::Initialize_UI()
{
	if(use_ui){
		ImguiGlut::Initialize_Imgui();
		////Initialize window properties
		ImGui::GetStyle().WindowRounding=0.0f;
		ImGui::GetStyle().FrameRounding=0.0f;
		ImGui::GetStyle().GrabRounding=0.0f;
		ImGui::GetStyle().ScrollbarRounding=0.0f;
		ImGui::GetStyle().FramePadding=ImVec2(0.0f,0.0f);
		ImGui::GetStyle().Alpha=.5f;
		ImGui::GetStyle().WindowTitleAlign=ImVec2(.5f,.5f);

		////Initialize command window and set it invisible
		ui_command.reset(new OpenGLUICommand());
		ui_command->opengl_viewer=this;
		ui_command->visible=false;
		Add_Object(ui_command.get());
	}
}

void OpenGLViewer::Finish()
{
	if(use_ui)ImguiGlut::Shut_Down_Imgui();
}

////sentences are separated by ';'
void Get_Sentences(const std::string& cmd,Array<std::string>& sentences)
{
	std::stringstream ss(cmd);std::string s;
	while(std::getline(ss,s,';')){sentences.push_back(s);}
}

////tokens are separated by ' '
void Get_Tokens(const std::string& cmd,Array<std::string>& tokens)
{
	std::stringstream ss(cmd);std::string s;
	while(std::getline(ss,s,' ')){if(s!=" "&&s!="")tokens.push_back(s);}
}

void OpenGLViewer::Toggle_Command(const std::string cmd)
{
	Array<std::string> sentences;Get_Sentences(cmd,sentences);
	for(int i=0;i<(int)sentences.size();i++){Toggle_Single_Command(sentences[i]);}
}

void OpenGLViewer::Toggle_Single_Command(const std::string cmd)
{
	//std::cout<<"cmd: "<<cmd<<std::endl;
	Array<std::string> tokens;
	Get_Tokens(cmd,tokens);
	if(tokens.empty())return;
	bool valid_cmd=false;

	//std::cout<<"token_size: "<<tokens.size()<<std::endl;
	//for(int i=0;i<tokens.size();i++)std::cout<<tokens[i]<<", ";std::cout<<std::endl;

	////slice view
	////Command: slice [string:grid_name] [int:axis] [int:start]
	if(tokens[0]=="slice"&&tokens.size()==4){
		std::string name=tokens[1];
		int axis=std::stoi(tokens[2]);
		int start=std::stoi(tokens[3]);
		std::cout<<"Cmd [slice]: name: "<<name<<", axis: "<<axis<<", start: "<<start<<std::endl;
		auto* obj=Get_Object_By_Name(name);
		if(obj==nullptr){std::cerr<<"Invalid cmd [slice]: cannot find object by name "<<name<<std::endl;return;}
		OpenGLGrid* grid_obj=dynamic_cast<OpenGLGrid*>(obj);
		if(grid_obj==nullptr){std::cerr<<"Invalid cmd [slice]: invalid grid object "<<name<<std::endl;return;}
		std::cout<<"Action: update grid: "<<grid_obj->name<<std::endl;
		Set_Subgrid_Slice(grid_obj,axis,start);
		grid_obj->Set_Data_Refreshed();

		for(auto& obj:opengl_window->object_list){
			auto opengl_object=dynamic_cast<OpenGLObject*>(obj.get());
			std::string obj_name=opengl_object->name;
			if(obj_name.find(name+" ")!=std::string::npos){
				OpenGLGridField<real>* grid_field_obj=dynamic_cast<OpenGLGridField<real>* >(opengl_object);
				if(grid_field_obj!=nullptr){
					std::cout<<"Action: update grid scalar field: "<<obj_name<<std::endl;
					Set_Subgrid_Slice(grid_field_obj,axis,start);
					grid_field_obj->Set_Data_Refreshed();
					continue;}

				OpenGLGridVectors<real>* grid_vectors_obj=dynamic_cast<OpenGLGridVectors<real>* >(opengl_object);
				if(grid_vectors_obj!=nullptr){
					std::cout<<"Action: update grid vector field: "<<obj_name<<std::endl;
					Set_Subgrid_Slice(grid_vectors_obj,axis,start);
					grid_vectors_obj->Set_Data_Refreshed();
					continue;}}}
		valid_cmd=true;}

	////slice view off
	////Command: slice [string:grid_name] off
	if(tokens.size()==3&&tokens[0]=="slice"&&tokens[2]=="off"){
		std::string name=tokens[1];
		std::cout<<"Cmd [slice]: name: "<<name<<" off"<<std::endl;
		auto* obj=Get_Object_By_Name(name);
		if(obj==nullptr){std::cerr<<"Invalid cmd [slice]: cannot find object by name "<<name<<std::endl;return;}
		OpenGLGrid* grid_obj=dynamic_cast<OpenGLGrid*>(obj);
		if(grid_obj==nullptr){std::cerr<<"Invalid cmd [slice]: invalid grid object "<<name<<std::endl;return;}
		std::cout<<"Action: update grid: "<<grid_obj->name<<std::endl;
		grid_obj->Use_Subgrid_View(false);
		grid_obj->Set_Data_Refreshed();

		for(auto& obj:opengl_window->object_list){
			auto opengl_object=dynamic_cast<OpenGLObject*>(obj.get());
			std::string obj_name=opengl_object->name;
			if(obj_name.find(name+" ")!=std::string::npos){
				OpenGLGridField<real>* grid_field_obj=dynamic_cast<OpenGLGridField<real>* >(opengl_object);
				if(grid_field_obj!=nullptr){
					std::cout<<"Action: update grid scalar field: "<<obj_name<<std::endl;
					grid_field_obj->Use_Subgrid_View(false);
					grid_field_obj->Set_Data_Refreshed();
					continue;}

				OpenGLGridVectors<real>* grid_vectors_obj=dynamic_cast<OpenGLGridVectors<real>* >(opengl_object);
				if(grid_vectors_obj!=nullptr){
					std::cout<<"Action: update grid vector field: "<<obj_name<<std::endl;
					grid_vectors_obj->Use_Subgrid_View(false);
					grid_vectors_obj->Set_Data_Refreshed();
					continue;}}}
		valid_cmd=true;}

	////goto frame
	////Command: goto [int:frame_idx]
	if(tokens.size()==2&&tokens[0]=="goto"){
		int target_frame=std::stoi(tokens[1]);
		frame=target_frame;
		valid_cmd=true;}

	////offscreen
	////Command: offscreen on/off
	if(tokens.size()==2&&tokens[0]=="offscreen"){
		if(tokens[1]=="on")opengl_window->display_offscreen=true;
		else if(tokens[1]=="off")opengl_window->display_offscreen=false;
		valid_cmd=true;}

	////set frame interval
	////Command: every [int:frame_interval]
	if(tokens.size()==2&&tokens[0]=="every"){
		int every_frame=std::stoi(tokens[1]);
		every=every_frame;
		valid_cmd=true;}

	////save camera
	if(tokens.size()==2&&tokens[0]=="save"&&tokens[1]=="camera"){
		//std::string file_name=output_dir+"/camera.txt";
		std::string file_name="camera.txt";
		opengl_window->Save_Camera(file_name);}

	////load camera
	if(tokens.size()==2&&tokens[0]=="load"&&tokens[1]=="camera"){
		//std::string file_name=output_dir+"/camera.txt";
		std::string file_name="camera.txt";
		opengl_window->Load_Camera(file_name);}

	////text off
	if(tokens.size()==2&&tokens[0]=="text"){
		if(tokens[1]=="on")opengl_window->display_text=true;
		else if(tokens[1]=="off")opengl_window->display_text=false;}

	////update frame after executing command
	if(valid_cmd)Update_Frame();
}

//////////////////////////////////////////////////////////////////////////
////Set viewer and object properties

void OpenGLViewer::Set_Offscreen_Output_Dir(const std::string offscreen_output)
{
	opengl_window->offscreen_output_dir=offscreen_output;
}

//////////////////////////////////////////////////////////////////////////
////Add objects

template<class T_FIELD> T_FIELD* Add_Grid_Field(OpenGLViewer& viewer,const std::string& grid_name,const Array<OpenGLData> data=Array<OpenGLData>())
{
	T_FIELD* opengl_grid_field=nullptr;
	Array<OpenGLData> readable_data;
	for(auto i=0;i<data.size();i++)
		if(OpenGLObject::Object_File_Exists(viewer.output_dir,viewer.first_frame,data[i].name)){readable_data.push_back(data[i]);}
	if(readable_data.size()>0&&viewer.Initialize_From_File(viewer.output_dir,opengl_grid_field,grid_name,readable_data,viewer.first_frame)){
		opengl_grid_field->Refresh(viewer.first_frame);
		viewer.Bind_Callback_Keys(readable_data,opengl_grid_field->callbacks,0);
		viewer.opengl_window->Add_Object(opengl_grid_field);
		for(int i=0;i<readable_data.size();i++)
			opengl_grid_field->name+=(" "+readable_data[i].name);
		if(viewer.verbose){std::cout<<"Add opengl grid field";AuxFunc::Seperation_Line();}}
	return opengl_grid_field;
}

//OpenGLVolume* OpenGLViewer::Add_Grid_Volume(const std::string& grid_name,const Array<OpenGLData> data)
//{
//	auto opengl_vol=Add_Object<OpenGLVolume>();
//
//	Array<OpenGLData> readable_data;
//	for(auto i=0;i<data.size();i++)
//		if(OpenGLObject::Object_File_Exists(output_dir,first_frame,data[i].name)){readable_data.push_back(data[i]);}
//	if(readable_data.size()>0&&Initialize_From_File(output_dir,opengl_vol,grid_name,readable_data,first_frame)){
//		opengl_vol->Refresh(first_frame);
//		Bind_Callback_Keys(readable_data,opengl_vol->callbacks,0);
//		opengl_window->Add_Object(opengl_vol);
//		if(verbose){std::cout<<"Add opengl grid volume";AuxFunc::Seperation_Line();}}
//	return opengl_vol;
//}

template<class T_MESH,class T_FIELD> T_FIELD* Add_Mesh_Field(OpenGLViewer& viewer,T_MESH* opengl_mesh,const Array<OpenGLData> data=Array<OpenGLData>())
{
	T_FIELD* opengl_field=nullptr;
	Array<OpenGLData> readable_data;
	for(auto i=0;i<data.size();i++)
		if(OpenGLObject::Object_File_Exists(viewer.output_dir,viewer.first_frame,data[i].name)){readable_data.push_back(data[i]);}
	if(opengl_mesh!=nullptr&&readable_data.size()>0){
		opengl_field=new T_FIELD(opengl_mesh);
		opengl_field->output_dir=viewer.output_dir;
		opengl_field->data=readable_data;
		opengl_field->Refresh(viewer.first_frame);
		viewer.Bind_Callback_Keys(readable_data,opengl_field->callbacks,0);
		viewer.opengl_window->Add_Object(opengl_field);
		if(viewer.verbose){std::cout<<"Add OpenGL mesh field";AuxFunc::Seperation_Line();}}
	return opengl_field;
}

OpenGLGridField<real>* OpenGLViewer::Add_Grid_Scalar_Field(const std::string& grid_name,const Array<OpenGLData> data/*=Array<OpenGLData>()*/)
{return Add_Grid_Field<OpenGLGridField<real> >(*this,grid_name,data);}

OpenGLGridVectors<real>* OpenGLViewer::Add_Grid_Vector_Field(const std::string& grid_name,const Array<OpenGLData> data/*=Array<OpenGLData>()*/)
{return Add_Grid_Field<OpenGLGridVectors<real> >(*this,grid_name,data);}

OpenGLGridTensors<real>* OpenGLViewer::Add_Grid_Tensor_Field(const std::string& grid_name,const Array<OpenGLData> data/*=Array<OpenGLData>()*/)
{return Add_Grid_Field<OpenGLGridTensors<real> >(*this,grid_name,data);}

OpenGLGridHeightField<real>* OpenGLViewer::Add_Grid_Height_Field(const std::string& grid_name,const Array<OpenGLData> data/*=Array<OpenGLData>()*/)
{return Add_Grid_Field<OpenGLGridHeightField<real> >(*this,grid_name,data);}

OpenGLVolume<real>* OpenGLViewer::Add_Grid_Volume(const std::string& grid_name,const Array<OpenGLData> data)
{return Add_Grid_Field<OpenGLVolume<real> >(*this,grid_name,data);}

OpenGLMeshField<real,TriangleMesh<3> >* OpenGLViewer::Add_Mesh_Scalar_Field(OpenGLTriangleMesh* opengl_mesh,const Array<OpenGLData> data)
{return Add_Mesh_Field<OpenGLTriangleMesh,OpenGLMeshField<real,TriangleMesh<3> > >(*this,opengl_mesh,data);}

OpenGLMeshField<real,TetrahedronMesh<3> >* OpenGLViewer::Add_Mesh_Scalar_Field(OpenGLTetrahedronMesh* opengl_mesh,const Array<OpenGLData> data)
{return Add_Mesh_Field<OpenGLTetrahedronMesh,OpenGLMeshField<real,TetrahedronMesh<3> > >(*this,opengl_mesh,data);}

OpenGLMeshVectors<real,TriangleMesh<3> >* OpenGLViewer::Add_Mesh_Vector_Field(OpenGLTriangleMesh* opengl_mesh,const Array<OpenGLData> data)
{return Add_Mesh_Field<OpenGLTriangleMesh,OpenGLMeshVectors<real,TriangleMesh<3> > >(*this,opengl_mesh,data);}

OpenGLMeshVectors<real,TetrahedronMesh<3> >* OpenGLViewer::Add_Mesh_Vector_Field(OpenGLTetrahedronMesh* opengl_mesh,const Array<OpenGLData> data)
{return Add_Mesh_Field<OpenGLTetrahedronMesh,OpenGLMeshVectors<real,TetrahedronMesh<3> > >(*this,opengl_mesh,data);}

OpenGLMeshTensors<real,TriangleMesh<3> >* OpenGLViewer::Add_Mesh_Tensor_Field(OpenGLTriangleMesh* opengl_mesh,const Array<OpenGLData> data)
{return Add_Mesh_Field<OpenGLTriangleMesh,OpenGLMeshTensors<real,TriangleMesh<3> > >(*this,opengl_mesh,data);}

OpenGLMeshTensors<real,TetrahedronMesh<3> >* OpenGLViewer::Add_Mesh_Tensor_Field(OpenGLTetrahedronMesh* opengl_mesh,const Array<OpenGLData> data)
{return Add_Mesh_Field<OpenGLTetrahedronMesh,OpenGLMeshTensors<real,TetrahedronMesh<3> > >(*this,opengl_mesh,data);}

OpenGLText* OpenGLViewer::Add_Interactive_Text_2D(const std::string& text,const Vector2& pos)
{
	OpenGLText* opengl_text=Add_Interactive_Object<OpenGLText>();
	opengl_text->texts.push_back(text);
	opengl_text->offsets.push_back(Vector2(0,0));
	auto camera=OpenGLUbos::Get_Camera_Ubo();
	opengl_text->Set_Data_Pointers(&camera->object.ortho);
	opengl_text->Set_Pos(pos);
	return opengl_text;
}

OpenGLText3D* OpenGLViewer::Add_Interactive_Text_3D(const std::string& text,const Vector3& pos)
{
	OpenGLText3D* opengl_text=Add_Interactive_Object<OpenGLText3D>();
	opengl_text->opengl_window=opengl_window.get();
	opengl_text->texts.push_back(text);
	opengl_text->offsets.push_back(Vector2(0,0));
	auto camera=OpenGLUbos::Get_Camera_Ubo();
	opengl_text->Set_Data_Pointers(&camera->object.ortho);
	opengl_text->Set_Pos(pos);
	return opengl_text;
}

OpenGLTextArray3D* OpenGLViewer::Add_Interactive_Text_Array_3D(const Array<std::string>& text,const Array<Vector3>& pos)
{
	OpenGLTextArray3D* opengl_text=Add_Interactive_Object<OpenGLTextArray3D>();
	opengl_text->opengl_window=opengl_window.get();
	opengl_text->texts=text;
	for(int i=0;i<text.size();i++)opengl_text->offsets.push_back(Vector2(0,0));
	opengl_text->pos_3d_array=pos;
	auto camera=OpenGLUbos::Get_Camera_Ubo();
	opengl_text->Set_Data_Pointers(&camera->object.ortho);
	return opengl_text;
}

OpenGLTriangleMesh* OpenGLViewer::Add_Interactive_Triangle_Mesh()
{
	return Add_Interactive_Object<OpenGLTriangleMesh>();
}

OpenGLParticles<Particles<3> >* OpenGLViewer::Add_Interactive_Particles()
{
	return Add_Interactive_Object<OpenGLParticles<Particles<3> > >();
}

OpenGLSphere* OpenGLViewer::Add_Interactive_Sphere(const Vector3& center/*=Vector3::Zero()*/,const real r/*=(real)1*/)
{
	auto opengl_sphere=Add_Interactive_Object<OpenGLSphere>();
	opengl_sphere->pos=center;
	opengl_sphere->radius=r;
	opengl_sphere->Initialize();
	return opengl_sphere;
}

OpenGLArrow* OpenGLViewer::Add_Interactive_Arrow(const Vector3& start/*=Vector3::Zero()*/, const Vector3& end/*=Vector3::Unit(1)*/)
{
	auto opengl_arrow=Add_Interactive_Object<OpenGLArrow>();
	opengl_arrow->start=start;
	opengl_arrow->end=end;
	opengl_arrow->Initialize();
	return opengl_arrow;
}

OpenGLObject* OpenGLViewer::Get_Object_By_Name(const std::string& name)
{
	for(auto& obj_ptr:opengl_window->object_list){
		OpenGLObject* obj=obj_ptr.get();
		if(obj->name==name)return obj;}
	return nullptr;
}

OpenGLBackground* OpenGLViewer::Get_Object_Background()
{
	OpenGLBackground* opengl_background=Get_Object_By_Name_And_Type<OpenGLBackground>("background");
	return opengl_background;
}

void OpenGLViewer::Set_Background_Color(const OpenGLColor& c1,const OpenGLColor& c2)
{
	auto opengl_background=Get_Object_Background();
	if(opengl_background){
		opengl_background->mix_colors[0]=c1;
		opengl_background->mix_colors[1]=c2;
		opengl_background->Set_Data_Refreshed();}
}

OpenGLTriangleMesh* OpenGLViewer::Get_Object_Triangle_Mesh_By_Name(const std::string& name)
{
	OpenGLTriangleMesh* opengl_triangle_mesh=Get_Object_By_Name_And_Type<OpenGLTriangleMesh>(name);
	return opengl_triangle_mesh;
}

//////////////////////////////////////////////////////////////////////////
////Basic callbacks

void OpenGLViewer::Print_Keyboard_Callbacks()
{
	AuxFunc::Seperation_Line();
	std::cout<<"Keyboard callback functions:\n";
	for(auto iter=key_data_hashtable.begin();iter!=key_data_hashtable.end();iter++){
		std::cout<<"Key: "<<iter->first<<", func: "<<iter->second<<std::endl;}
	AuxFunc::Seperation_Line();
}

void OpenGLViewer::Object_Info_Debug()
{
	AuxFunc::Seperation_Line();
	std::cout<<"OpenGL objects:\n";
	for (auto& obj : opengl_window->object_list) {
		obj->Output_Debug_Info(std::cout);
		std::cout << std::endl;
	}
	AuxFunc::Seperation_Line();
	opengl_window->Toggle_Color_Encoding();
}

void OpenGLViewer::Toggle_Increase_Scale()
{
	for(auto& obj:opengl_window->object_list){
		auto opengl_object=dynamic_cast<OpenGLObject*>(obj.get());
		if(opengl_object->visible)opengl_object->Toggle_Increase_Scale();}
}

void OpenGLViewer::Toggle_Decrease_Scale()
{
	for(auto& obj:opengl_window->object_list){
		auto opengl_object=dynamic_cast<OpenGLObject*>(obj.get());
		if(opengl_object->visible)opengl_object->Toggle_Decrease_Scale();}
}

void OpenGLViewer::Toggle_Increase_Feature()
{
	for (auto& obj : opengl_window->object_list) {
		auto opengl_object = dynamic_cast<OpenGLObject*>(obj.get());
		if (opengl_object->visible)opengl_object->Toggle_Increase_Feature();
	}
}

void OpenGLViewer::Toggle_Decrease_Feature()
{
	for (auto& obj : opengl_window->object_list) {
		auto opengl_object = dynamic_cast<OpenGLObject*>(obj.get());
		if (opengl_object->visible)opengl_object->Toggle_Decrease_Feature();
	}
}

void OpenGLViewer::Toggle_Normalize_Data()
{
	for(auto& obj:opengl_window->object_list){
		auto opengl_object=dynamic_cast<OpenGLObject*>(obj.get());
		if(opengl_object->visible)opengl_object->Normalize_Data();}
}

void OpenGLViewer::Toggle_Draw_Dis()
{
	for(auto& obj:opengl_window->object_list){
		auto opengl_object=dynamic_cast<OpenGLObject*>(obj.get());
		if(opengl_object->visible)opengl_object->Toggle_Draw_Dis();}
}

void OpenGLViewer::Toggle_Next_Frame()
{
	frame=frame+every;
	if(last_frame!=-1&&frame>last_frame){
		if(frame>last_frame){	////refresh last frame
			frame=last_frame-(last_frame%every);
			std::string file_name=output_dir+"/0/last_frame.txt";
			if(File::File_Exists(file_name)){
				int new_last_frame=last_frame;
				File::Read_Text_From_File(file_name,new_last_frame);
				if(new_last_frame>last_frame)last_frame=new_last_frame-(new_last_frame%every);}
			return;}}
	Update_Frame();
}

void OpenGLViewer::Toggle_Prev_Frame()
{
	frame=frame-every;
	if(frame<0)frame=0;
	else Update_Frame();
}

void OpenGLViewer::Toggle_First_Frame()
{
	frame=first_frame;
	Update_Frame();
}

void OpenGLViewer::Toggle_Play()
{
	play=!play;
	opengl_window->Set_Timer_Callback(play?&Toggle_Next_Frame_Func:nullptr);
}

void OpenGLViewer::Initialize_Common_Callback_Keys()
{
	Bind_Callback_Key('>',&Toggle_Increase_Scale_Func,"scale+");
	Bind_Callback_Key('<',&Toggle_Decrease_Scale_Func,"scale-");
	Bind_Callback_Key('=', &Toggle_Increase_Feature_Func, "feature+");
	Bind_Callback_Key('-', &Toggle_Decrease_Feature_Func, "feature-");
	Bind_Callback_Key('K',&Print_Keyboard_Callbacks_Func,"print binded keys");
	Bind_Callback_Key('O',&Object_Info_Debug_Func,"print object names and debug");
	Bind_Callback_Key(']',&Toggle_Next_Frame_Func,"next frame");
	Bind_Callback_Key('[',&Toggle_Prev_Frame_Func,"prev frame");
	Bind_Callback_Key('r',&Toggle_First_Frame_Func,"first frame");
	Bind_Callback_Key('p',&Toggle_Play_Func,"play");
	Bind_Callback_Key('q',&opengl_window->Quit_Func,"quit");
	Bind_Callback_Key('w',&opengl_window->Toggle_Offscreen_Func,"offscreen rendering");
}

void OpenGLViewer::Bind_Callback_Key(const uchar key,std::function<void(void)>* callback,const std::string& discription)
{
	opengl_window->Set_Keyboard_Callback(std::string(1,key),callback);
	key_data_hashtable.insert(std::make_pair(key,discription));
}

void OpenGLViewer::Bind_Callback_Keys(const Array<OpenGLData>& data,Array<std::function<void(void)>*> data_idx_callbacks,int start_idx/*=0*/)
{
	for(size_type i=0;i<data.size();i++){if(data[i].key.size()>0){
		Bind_Callback_Key(data[i].key[0],data_idx_callbacks[i+start_idx],data[i].name);}}
}
