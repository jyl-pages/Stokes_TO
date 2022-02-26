//////////////////////////////////////////////////////////////////////////
// Opengl world
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "OpenGLWindow.h"
#include <iostream>
#include <GL/freeglut.h>
#include "glm.hpp"
#include "gtc/matrix_transform.hpp"
#include "gtc/type_ptr.hpp"
#include "StbImage.h"
#include "File.h"
#include "GeometryPrimitives.h"
#include "OpenGLColor.h"
#include "OpenGLArcball.h"
#include "OpenGLObject.h"
#include "OpenGLFbos.h"
#include "OpenGLUbos.h"
#include "OpenGLShaderLibrary.h"
#include "OpenGLViewer.h"
//#include "OpenGLCamera.h"

using namespace OpenGLUbos;
using namespace OpenGLFbos;

OpenGLWindow* OpenGLWindow::instance=nullptr;
const unsigned int time_per_frame=20;

////////////////////////////////////////////////////////////////////////////////////////////////////
OpenGLWindow::OpenGLWindow()
{
	instance=this;

	arcball=std::make_shared<OpenGLArcball>();
	arcball->window=this;
	arcball_matrix=arcball->Value();
	rotation_matrix=arcball->Value();
}

void OpenGLWindow::Init()
{
	Initialize_Window(); // Initialize window, and set a set of callback functions
	Initialize_OpenGL();
}

void OpenGLWindow::Run()
{
	glutMainLoop();
}

void OpenGLWindow::Initialize_Window()
{
	int argc=1;char* argv[1];argv[0]=(char*)(window_title.c_str());

	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGBA|GLUT_DEPTH|GLUT_ALPHA);
	glutInitWindowSize(win_w,win_h);
	window_id=glutCreateWindow(window_title.c_str());

	glutSetWindow(window_id);
	glutIdleFunc(Idle_Func_Glut); // It's generally nothing but calls glutPostRedisplay() at the end of every main loop
	glutTimerFunc(time_per_frame,Timer_Func_Glut,1); // It calls OpenGLViewer::Toggle_Next_Frame() if play is on.
	glutDisplayFunc(Display_Func_Glut); // Called at the beginning of every main loop. It's basically OpenGLWindow::Display()
	glutReshapeFunc(Reshape_Func_Glut);
	glutMouseFunc(Mouse_Func_Glut);
	glutMotionFunc(Motion_Func_Glut);
	glutKeyboardFunc(Keyboard_Func_Glut);
	glutKeyboardUpFunc(Keyboard_Up_Func_Glut);
	glutSpecialFunc(Keyboard_Special_Func_Glut);
	glutSpecialUpFunc(Keyboard_Special_Up_Func_Glut);
}

void OpenGLWindow::Initialize_OpenGL()
{
	// Get and output OpenGL version info
	GLint major_version;glGetIntegerv(GL_MAJOR_VERSION,&major_version);
	GLint minor_version;glGetIntegerv(GL_MINOR_VERSION,&minor_version);
	std::cout<<"Opengl major version: "<<major_version<<", minor version: "<<minor_version<<std::endl;

	glutInitContextVersion(4,3);
	glutInitContextProfile(GLUT_CORE_PROFILE);
	if(glewInit()){std::cerr<<"Error: [OpenGLWindow] Cannot initialize glew"<<std::endl;return;}

	glEnable(GL_DEPTH_TEST); // Do use depth buffer
//	glPixelStorei(GL_PACK_ALIGNMENT,1);
//	glPixelStorei(GL_UNPACK_ALIGNMENT,1);
	glFrontFace(GL_CCW); // Counter-clockwise as front side

	Initialize_Camera();
	Initialize_Ubos(); // It's in OpenGLUbos.h
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void OpenGLWindow::Display()
{
	Update_Camera(); // Update camera parameters and bind it to uniform buffer
	Preprocess();
	Clear_Buffers();
	for(auto& obj:object_list){obj->Display();}

	if(display_text)Display_Text();
	if(display_offscreen)Display_Offscreen();

	GLenum gl_error=glGetError();
	if(gl_error!=GL_NO_ERROR){std::cerr<<"Error: [OpenGLWindow] "<<gluErrorString(gl_error)<<std::endl;}
}

void OpenGLWindow::Preprocess()
{
	bool use_preprocess=false;for(auto& obj:object_list){
		OpenGLObject* o=dynamic_cast<OpenGLObject*>(obj.get());
		if(o->use_preprocess){use_preprocess=true;break;}}if(!use_preprocess)return;

	bool use_depth_fbo=false;for(auto& obj:object_list){
		OpenGLObject* o=dynamic_cast<OpenGLObject*>(obj.get());
		if(o->use_depth_fbo){use_depth_fbo=true;break;}}
	if(use_depth_fbo){
		GLuint depth_w=1024,depth_h=1024;
		auto fbo=Get_Depth_Fbo("depth");
		fbo->Resize(depth_w,depth_h);fbo->Clear();}

	for(auto& obj:object_list){
		OpenGLObject* o=dynamic_cast<OpenGLObject*>(obj.get());if(!o->use_preprocess)continue;
		o->Preprocess();}

	//if(use_depth_fbo){
	//	auto fbo=Get_Depth_Fbo("depth");
	//	fbo->Write_To_File("depth_map");}
}

void OpenGLWindow::Update_Data_To_Render()
{
	for(auto& obj:object_list){obj->Update_Data_To_Render();}
}

void OpenGLWindow::Redisplay()
{
	glutPostRedisplay();
}

void OpenGLWindow::Display_Offscreen()
{
    ////render to image
    if((!display_offscreen_interactive&&frame_offscreen!=frame_offscreen_rendered)||display_offscreen_interactive){
        int w=win_w;int h=win_h;
        int num_pixel=w*h;int num_comp=3;
        GLubyte* pixels=new GLubyte[num_comp*num_pixel];
        GLubyte* pixels_flipped_y=new GLubyte[num_comp*num_pixel];
        glReadPixels(0,0,w,h,GL_RGB,GL_UNSIGNED_BYTE,pixels);
        for(int i=0;i<h;i++){int offset=w*num_comp*(h-i-1);
            std::memcpy(pixels_flipped_y+offset,pixels+w*num_comp*i,w*num_comp);}
        int wrt_frame=display_offscreen_interactive?frame_offscreen_rendered++:frame_offscreen;
        if(!File::Directory_Exists(offscreen_output_dir.c_str()))File::Create_Directory(offscreen_output_dir);
        std::stringstream ss;ss<<offscreen_output_dir<<"/"<<std::setfill('0')<<std::setw(4)<<wrt_frame<<".png";
        std::cout<<"Offscreen render to image "<<ss.str()<<std::endl;
        int rst=Stb::Write_Png(ss.str().c_str(),w,h,num_comp,pixels_flipped_y,0);
        delete pixels;delete pixels_flipped_y;
        if(!display_offscreen_interactive)frame_offscreen_rendered=frame_offscreen;}
}

void OpenGLWindow::Display_Text()
{
	if(texts.empty())return;

    auto camera=Get_Camera_Ubo();
    glm::mat4& ortho=camera->object.ortho;

    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadMatrixf(glm::value_ptr(ortho));
    glColor3f(.5f,.5f,.5f);

    Vector2f step={10.f,15.f};
	int i=0;for(auto iter:texts){const std::string& text=iter.second;
		for(auto j=0;j<text.size();j++){
			Vector2f pos=Vector2f(step[0]*(float)j,(float)win_h-step[1]*((float)(i)+1.f));
			glRasterPos2f(pos[0],pos[1]);
			glutBitmapCharacter(GLUT_BITMAP_9_BY_15,text[j]);}i++;}

    glPopAttrib();
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void OpenGLWindow::Add_Object(OpenGLObject *object)
{
    object_list.push_back(std::unique_ptr<OpenGLObject>(object));
}

void OpenGLWindow::Add_Object(std::unique_ptr<OpenGLObject>& object)
{
    object_list.push_back(std::move(object));
}

void OpenGLWindow::Clear_Objects()
{
    object_list.clear();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void OpenGLWindow::Initialize_Camera()
{
    camera_target=Vector3f::Zero();
    camera_distance=10.f;
	// Updates nearclip and farclip
    Update_Clip_Planes();
}

void OpenGLWindow::Update_Camera()
{
    auto camera=Get_Camera_Ubo(); // It'a pointer to a OpenGLUboInstance<Camera>
    glm::mat4& proj=camera->object.projection;
    proj=glm::perspective(glm::radians(fovy),(float)win_w/(float)win_h,nearclip,farclip);

    glm::mat4& view=camera->object.view;
    if(false/*splined_camera.initialized*/){
        /*Vector3f pos=splined_camera.Get_Pos().cast<float>();
        Vector3f dir=splined_camera.Get_Dir().cast<float>();
        view=glm::lookAt(glm::vec3(pos[0],pos[1],pos[2]),glm::vec3(dir[0],dir[1],dir[2]),glm::vec3(.0f,1.f,.0f));*/}
    else{
        view=glm::translate(glm::mat4(),glm::vec3(0.f,0.f,(float)-camera_distance))*glm::make_mat4x4(arcball_matrix.data())*glm::make_mat4x4(rotation_matrix.data());
        view=glm::translate(view,glm::vec3(-camera_target[0],-camera_target[1],-camera_target[2]));}

    glm::mat4& pvm=camera->object.pvm;
    pvm=proj*view;	////assuming model matrix is identity
    glm::mat4& ortho=camera->object.ortho;
    ortho=glm::ortho(0.f,(GLfloat)win_w,0.f,(GLfloat)win_h);
    glm::vec4& position=camera->object.position;
    glm::mat4 inv_view=glm::inverse(view);
    position=glm::vec4(inv_view[3][0],inv_view[3][1],inv_view[3][2],1.f);
    camera->Set_Block_Attributes();
}

void OpenGLWindow::Save_Camera(const std::string camera_file_name)
{
	////save camera_distance, arcball_matrix,rotation_matrix,camera target, 36 doubles
	Array<real> camera_data;
	camera_data.push_back(camera_distance);
	for(int i=0;i<16;i++)camera_data.push_back(arcball_matrix.data()[i]);
	for(int i=0;i<16;i++)camera_data.push_back(rotation_matrix.data()[i]);
	for(int i=0;i<3;i++)camera_data.push_back(camera_target[i]);
	File::Write_Text_Array_To_File(camera_file_name,camera_data,36);
	std::cout<<"Write camera data to file "<<camera_file_name<<std::endl;
}

void OpenGLWindow::Load_Camera(const std::string camera_file_name)
{
	if(!File::File_Exists(camera_file_name)){
		std::cerr<<"[Error] Camera file does not exist: "<<camera_file_name<<std::endl;
		return;}

	////read camera_distance, arcball_matrix,rotation_matrix,camera target
	Array<real> camera_data;int size=36;camera_data.resize(size);
	File::Read_Text_Array_From_File(camera_file_name,camera_data,36);
	//std::cout<<"camera_data: "<<camera_data.size()<<": ";
	//for(int i=0;i<camera_data.size();i++){std::cout<<camera_data[i]<<", ";}
	int p=0;
	camera_distance=camera_data[p++];
	for(int i=0;i<16;i++)arcball_matrix.data()[i]=camera_data[p++];
	for(int i=0;i<16;i++)rotation_matrix.data()[i]=camera_data[p++];
	for(int i=0;i<3;i++)camera_target.data()[i]=camera_data[p++];
	std::cout<<"Read camera data from file "<<camera_file_name<<std::endl;
	Update_Camera();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void OpenGLWindow::Idle_Func_Glut()
{
	instance->Idle_Func();
	glutPostRedisplay();
}

void OpenGLWindow::Timer_Func_Glut(int value)
{
	instance->Timer_Func();
	glutTimerFunc(time_per_frame,Timer_Func_Glut,1);
}

void OpenGLWindow::Display_Func_Glut()
{
	instance->Display();
	glutSwapBuffers();
}

void OpenGLWindow::Reshape_Func_Glut(int w,int h)
{
	instance->Reshape_Func(w,h);
}

void OpenGLWindow::Mouse_Func_Glut(int button,int state,int x,int y)
{
	instance->Mouse_Func(button,state,x,y);
}

void OpenGLWindow::Motion_Func_Glut(int x,int y)
{
	instance->Motion_Func(x,y);
}

void OpenGLWindow::Keyboard_Func_Glut(unsigned char key,int x,int y)
{
	instance->Keyboard_Func(key,x,y);
}

void OpenGLWindow::Keyboard_Up_Func_Glut(unsigned char key,int x,int y)
{
	instance->Keyboard_Up_Func(key,x,y);
}

void OpenGLWindow::Keyboard_Special_Func_Glut(int key,int x,int y)
{
	instance->Keyboard_Special_Func(key,x,y);
}

void OpenGLWindow::Keyboard_Special_Up_Func_Glut(int key,int x,int y)
{
	instance->Keyboard_Special_Up_Func(key,x,y);
}

void OpenGLWindow::Idle_Func()
{
	if(idle_callback!=nullptr){(*idle_callback)();}
}

void OpenGLWindow::Timer_Func()
{
	if(timer_callback!=nullptr){(*timer_callback)();}
}

void OpenGLWindow::Reshape_Func(int w,int h)
{
    win_w=w;win_h=h;
    glViewport(0,0,(GLsizei)win_w,(GLsizei)win_h);
}

void OpenGLWindow::Mouse_Func(int button,int state,int x,int y)
{
	int bl_x=x;int bl_y=win_h-y;
	int left=(button==GLUT_LEFT_BUTTON?1:0)*(state==GLUT_DOWN?1:-1);
	int right=(button==GLUT_RIGHT_BUTTON?1:0)*(state==GLUT_DOWN?1:-1);
	int middle=(button==GLUT_MIDDLE_BUTTON?1:0)*(state==GLUT_DOWN?1:-1);

	bool is_focus_captured=false;
	for(auto& obj:object_list){
		is_focus_captured=obj->Mouse_Click(left,right,middle,bl_x,bl_y,Win_Width(),Win_Height());
		if(is_focus_captured){Redisplay();return;}}

	if(opengl_viewer!=nullptr){
		is_focus_captured=opengl_viewer->Mouse_Click(left,right,middle,bl_x,bl_y,Win_Width(),Win_Height());
		if(is_focus_captured){Redisplay();return;}}

	if (color_encoding_mode && state == GLUT_DOWN) {
		unsigned char data[4];
		glReadPixels(x, Win_Height() - 1 - y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, data);
		if (data[0] == 255 && data[1] == 255 && data[2] == 255) {
			Color_Code_Unselect();
		}
		else {
			int pickedID = data[0] + data[1] * 256 + data[2] * 256 * 256;
			Color_Code_Select(pickedID);
		}
	}

    switch(button){
        case GLUT_LEFT_BUTTON:
			if(use_2d_display)break;
            {Vector2f mouse_pos=Win_To_Norm_Coord(x,y);
                if(state==GLUT_UP){
                    if(mouse_state==MouseState::Rotation){
                        arcball->End_Drag(mouse_pos);
                        arcball_matrix=arcball->Value();
                        rotation_matrix=arcball_matrix*rotation_matrix;
                        arcball_matrix=Matrix4f::Identity();}
                    mouse_state=MouseState::None;}
                else if(state==GLUT_DOWN){
                    arcball->Begin_Drag(mouse_pos);
                    mouse_state=MouseState::Rotation;}
            }break;
        case GLUT_RIGHT_BUTTON:
            {if(state==GLUT_UP){mouse_state=MouseState::None;}
                else if(state==GLUT_DOWN){mouse_state=MouseState::Zoom;}
            }break;
        case GLUT_MIDDLE_BUTTON:
            {if(state==GLUT_UP){mouse_state=MouseState::None;}
                else if(state==GLUT_DOWN){mouse_state=MouseState::Motion;Update_Mouse_Drag_Target();}
            }break;}
    mouse_x=x;mouse_y=y;
}

void OpenGLWindow::Motion_Func(int x,int y)
{
	int bl_x=x;int bl_y=win_h-y;
	bool is_focus_captured=false;
	for(auto& obj:object_list){
		is_focus_captured=obj->Mouse_Drag(bl_x,bl_y,Win_Width(),Win_Height());
		if(is_focus_captured){Redisplay();return;}}

	if(opengl_viewer!=nullptr){
		is_focus_captured=opengl_viewer->Mouse_Drag(bl_x,bl_y,Win_Width(),Win_Height());
		if(is_focus_captured){Redisplay();return;}}

    switch(mouse_state){
        case MouseState::Rotation:
        {arcball->Update(Win_To_Norm_Coord(x,y));
            arcball_matrix=arcball->Value();}break;
        case MouseState::Zoom:
        {camera_distance*=pow(1.01f,-1.f*(float)(y-mouse_y));
            Update_Clip_Planes();}break;
        case MouseState::Motion:
        {float dx=(float)(mouse_x-x);float dy=(float)(y-mouse_y);
            camera_target+=dx*target_x_drag_vector+dy*target_y_drag_vector;}break;}
    mouse_x=x;mouse_y=y;
    Redisplay();
}

void OpenGLWindow::Keyboard_Func(unsigned char key,int x,int y)
{
	int bl_x=x;int bl_y=win_h-y;
	bool is_focus_captured=false;
	for(auto& obj:object_list){
		is_focus_captured=obj->Keyboard(key,bl_x,bl_y,true);
		if(is_focus_captured)break;}
	if(is_focus_captured){Redisplay();return;}

	////Check if a keyboard callback is available
    auto c=keyboard_callbacks.find(std::string(1,key));
    if(c!=keyboard_callbacks.end()){(*c->second)();Redisplay();}
}

void OpenGLWindow::Keyboard_Up_Func(unsigned char key,int x,int y)
{
	int bl_x=x;int bl_y=win_h-y;
	bool is_focus_captured=false;
	for(auto& obj:object_list){
		is_focus_captured=obj->Keyboard(key,bl_x,bl_y,false);
		if(is_focus_captured)break;}
	if(is_focus_captured){Redisplay();return;}
}

void OpenGLWindow::Keyboard_Special_Func(int key,int x,int y)
{
	int bl_x=x;int bl_y=win_h-y;
	bool is_focus_captured=false;
	for(auto& obj:object_list){
		is_focus_captured=obj->Keyboard_Special(key,bl_x,bl_y,true);
		if(is_focus_captured)break;}
	if(is_focus_captured){Redisplay();return;}
}

void OpenGLWindow::Keyboard_Special_Up_Func(int key,int x,int y)
{
	int bl_x=x;int bl_y=win_h-y;
	bool is_focus_captured=false;
	for(auto& obj:object_list){
		is_focus_captured=obj->Keyboard_Special(key,bl_x,bl_y,false);
		if(is_focus_captured)break;}
	if(is_focus_captured){Redisplay();return;}
}

void OpenGLWindow::Quit()
{
	exit(0);
}

void OpenGLWindow::Toggle_Offscreen()
{
	display_offscreen=!display_offscreen;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
Vector3f OpenGLWindow::Project(const Vector3f& pos)
{
    auto camera=Get_Camera_Ubo();
    glm::mat4 mv=camera->object.view;	////assuming model matrix is identity
    glm::mat4 proj=camera->object.projection;
    GLint viewport[4];glGetIntegerv(GL_VIEWPORT,viewport);
    glm::vec3 v(pos[0],pos[1],pos[2]);
    glm::vec4 vp;for(int i=0;i<4;i++)vp[i]=(GLfloat)viewport[i];
    glm::vec3 wpos=glm::project(v,mv,proj,vp);
    return Vector3f(wpos[0],wpos[1],wpos[2]);
}

Vector3f OpenGLWindow::Unproject(const Vector3f& win_pos)
{
    auto camera=Get_Camera_Ubo();
    glm::mat4 mv=camera->object.view;	////assuming model matrix is identity
    glm::mat4 proj=camera->object.projection;
    GLint viewport[4];glGetIntegerv(GL_VIEWPORT,viewport);
    glm::vec4 vp;for(int i=0;i<4;i++)vp[i]=(GLfloat)viewport[i];
    glm::vec3 world=glm::unProject(glm::vec3(win_pos[0],win_pos[1],win_pos[2]),mv,proj,vp);return Vector3f(world[0],world[1],world[2]);
}

GLfloat OpenGLWindow::Win_Depth(int bl_win_x,int bl_win_y)
{
	GLfloat depth;glReadPixels(bl_win_x,bl_win_y,1,1,GL_DEPTH_COMPONENT,GL_FLOAT,&depth);return depth;
}

Vector3f OpenGLWindow::Win_Coord_To_World_Coord(int bl_win_x,int bl_win_y)
{
	GLfloat depth=Win_Depth(bl_win_x,bl_win_y);
	//std::cout<<"depth: "<<bl_win_x<<", "<<bl_win_y<<": "<<depth<<std::endl;
	Vector3f win_pos=Vector3f((float)bl_win_x,(float)bl_win_y,depth);
	Vector3f pos=Unproject(win_pos);
	//std::cout<<"win: "<<win_pos.transpose()<<", depth: "<<depth<<", pos: "<<pos.transpose()<<std::endl;
	return pos;
}

// PhysBAM Function: This function converts mouse space pixel coordinates to the normalize coordinates the arcball expects
Vector2f OpenGLWindow::Win_To_Norm_Coord(int x,int y)
{
    if(win_w>=win_h){return Vector2f(((float)x/win_w-0.5f)*2.f*((float)win_w/(float)win_h),-((float)y/(float)win_h-0.5f)*2.f);}
    else{return Vector2f(((float)x/(float)win_w-0.5f)*2.f,-((float)y/(float)win_h-0.5f)*2.f*(win_h/win_w));}
}

void OpenGLWindow::Clear_Buffers()
{
    glClearColor(.7f,.7f,.7f,0.f);	////background color
    glClearDepth(1);
    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
}

void OpenGLWindow::Update_Clip_Planes()
{
	// Updates nearclip and farclip
    nearclip=.1f*camera_distance;
    farclip=10.f*camera_distance;
}

void OpenGLWindow::Update_Mouse_Drag_Target()
{
    Vector3f win_pos=Project(camera_target);
    Vector3f drag_target=Unproject(win_pos+Vector3f::Unit(0));
    target_x_drag_vector=drag_target-camera_target;
    drag_target=Unproject(win_pos+Vector3f::Unit(1));
    target_y_drag_vector=Vector3f(drag_target[0],drag_target[1],drag_target[2])-camera_target;
}

void OpenGLWindow::Update_Color_Code_Offset()
{
	int current_offset = 1;//0 is reserved for selecting nothing
	for (int i = 0; i < object_list.size(); i++) {
		auto& obj = object_list[i];
		obj->color_code_offset = current_offset;
		current_offset += obj->Color_Code_Size();
	}
	all_color_code_size = current_offset - 1;
}

void OpenGLWindow::Toggle_Color_Encoding()
{
	if (color_encoding_mode) {
		for (auto &obj : object_list) {
			obj->Recover_Shading_Mode();
		}
		color_encoding_mode = false;
		glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
	}
	else {
		for (auto& obj : object_list) {
			obj->Set_Shading_Mode(ShadingMode::ColorEncoding);
		}
		color_encoding_mode = true;
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	}
	Redisplay();
}

void OpenGLWindow::Color_Code_Select(int all_idx)
{
	if (all_color_code_size >= 256 * 256 * 256) {
		std::cout << "CAUTION: all " << all_color_code_size << " elements, possibly color encoding overflow\n";
	}
	for (int i = 0; i < object_list.size(); i++) {
		auto& obj = object_list[i];
		obj->Unselect_Element();
	}
	for (int i = 0; i < object_list.size(); i++) {
		auto& obj = object_list[i];
		int offset = obj->color_code_offset;
		if (offset <= all_idx && all_idx < offset + obj->Color_Code_Size()) {
			obj->Select_Element(all_idx - offset);
			break;
		}
	}
	Redisplay();
}

void OpenGLWindow::Color_Code_Unselect()
{
	for (auto& obj : object_list) {
		obj->Unselect_Element();
	}
	Redisplay();
}

void OpenGLWindow::Set_Keyboard_Callback(const std::string& key,std::function<void(void)>* callback)
{
    auto c=keyboard_callbacks.find(key);
    if(c==keyboard_callbacks.end()){
        keyboard_callbacks.insert(std::pair<std::string,std::function<void(void)>* >(key,callback));}
    else{keyboard_callbacks[key]=callback;}
}

void OpenGLWindow::Set_Idle_Callback(std::function<void(void)>* callback)
{
	idle_callback=callback;
}

void OpenGLWindow::Set_Timer_Callback(std::function<void(void)>* callback)
{
	timer_callback=callback;
}

GLuint Win_Width()
{
    if(OpenGLWindow::instance==nullptr)return -1;
    else return OpenGLWindow::instance->win_w;
}

GLuint Win_Height()
{
    if(OpenGLWindow::instance==nullptr)return -1;
    else return OpenGLWindow::instance->win_h;
}