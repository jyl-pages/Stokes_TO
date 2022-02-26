//////////////////////////////////////////////////////////////////////////
// Opengl world
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLWorld_h__
#define __OpenGLWorld_h__
#include <iomanip>
#include <map>
#include <functional>
#include <GL/glew.h>
#include "Common.h"
#include "OpenGLCommon.h"

////Forward declaration
class OpenGLObject;
class OpenGLArcball;
class OpenGLViewer;

class OpenGLWindow
{
public:
	static OpenGLWindow* instance;
	//// Glut window
	int window_id=0;
	std::string window_title="Simplex OpenGL Viewer";
	int win_w=1280,win_h=960;
	float fovy=30.f;

	//// Offscreen rendering
	bool display_offscreen=false;
	bool display_offscreen_interactive=false;
	int frame_offscreen=0;
	int frame_offscreen_rendered=-1;
	std::string offscreen_output_dir="offscreen_output";

	////Viewer
	std::shared_ptr<OpenGLViewer> opengl_viewer;
	////Objects
	Array<std::unique_ptr<OpenGLObject> > object_list;

	////Texts
	bool display_text=true;
	std::map<std::string,std::string> texts;

	////Camera
	Vector3f camera_target;
	float camera_distance=1.f;
	float nearclip,farclip;
	std::shared_ptr<OpenGLArcball> arcball;
	Matrix4f arcball_matrix;
	Matrix4f rotation_matrix;
	//OpenGLCamera splined_camera;
	
	////Interaction
	int mouse_x=0,mouse_y=0;	////store the old mouse pos
	enum class MouseState:int{None=0,Rotation,Zoom,Motion} mouse_state=MouseState::None;
	Vector3f target_x_drag_vector;
	Vector3f target_y_drag_vector;

	////Callbacks
	std::map<std::string,std::function<void(void)>* > keyboard_callbacks;
	std::function<void(void)>* idle_callback=nullptr; // Set in OpenGLWindow::Set_Idle_Callback()
	std::function<void(void)>* timer_callback=nullptr; // Set in OpenGLWindow::Set_Timer_Callback()

	////Dimension
	bool use_2d_display=false;

	////Color Encoding
	int all_color_code_size;
	bool color_encoding_mode=false;

public:
	OpenGLWindow();

	////Initialization and run
	void Init();
	void Run();
	void Initialize_Window();
	void Initialize_OpenGL();

	////Display
	void Display();
	void Preprocess();
	void Update_Data_To_Render();
	void Redisplay();
	void Display_Offscreen();

	////Objects
	void Add_Object(OpenGLObject* object);
	void Add_Object(std::unique_ptr<OpenGLObject>& object);
	void Clear_Objects();

	////Text
	void Display_Text();

	////Camera
	void Initialize_Camera();
	void Update_Camera();
	void Save_Camera(const std::string camera_file_name);
	void Load_Camera(const std::string camera_file_name);

	////Glut callbacks
	static void Idle_Func_Glut();
	static void Timer_Func_Glut(int value);
	static void Display_Func_Glut();
	static void Reshape_Func_Glut(int w,int h);
	static void Mouse_Func_Glut(int button,int state,int x,int y);
	static void Motion_Func_Glut(int x,int y);
	static void Keyboard_Func_Glut(unsigned char key,int x,int y);
	static void Keyboard_Up_Func_Glut(unsigned char key,int x,int y);
	static void Keyboard_Special_Func_Glut(int key,int x,int y);
	static void Keyboard_Special_Up_Func_Glut(int key,int x,int y);
	void Idle_Func();
	void Timer_Func();
	void Reshape_Func(int w,int h);
    void Mouse_Func(int button,int state,int x,int y);
    void Motion_Func(int x,int y);
	void Keyboard_Func(unsigned char key,int x,int y);
	void Keyboard_Up_Func(unsigned char key,int x,int y);
	void Keyboard_Special_Func(int key,int x,int y);
	void Keyboard_Special_Up_Func(int key,int x,int y);

	////Window-related callbacks	
	void Quit();
	Define_Function_Object(OpenGLWindow,Quit);
	void Toggle_Offscreen();
	Define_Function_Object(OpenGLWindow,Toggle_Offscreen);
	////Idle callback
	void Set_Idle_Callback(std::function<void(void)>* callback);
	////Timer callback
	void Set_Timer_Callback(std::function<void(void)>* callback); // It's only set in OpenGLViewer::Toggle_Play()

	////Keyboard callback
	void Set_Keyboard_Callback(const std::string& key,std::function<void(void)>* callback);

	////Helper functions
	Vector3f Project(const Vector3f& pos);
	Vector3f Unproject(const Vector3f& win_pos);
	GLfloat Win_Depth(int bl_win_x,int bl_win_y);
	Vector3f Win_Coord_To_World_Coord(int bl_win_x,int bl_win_y);	
	Vector2f Win_To_Norm_Coord(int x,int y);
	void Clear_Buffers();
	void Update_Clip_Planes();
	void Update_Mouse_Drag_Target();

	////Color encoding functions
	void Update_Color_Code_Offset();
	void Toggle_Color_Encoding();
	void Color_Code_Select(int all_idx);
	void Color_Code_Unselect();
};

////Global helper functions
GLuint Win_Width();
GLuint Win_Height();

#endif