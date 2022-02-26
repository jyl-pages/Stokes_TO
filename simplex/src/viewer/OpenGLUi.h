//////////////////////////////////////////////////////////////////////////
// Opengl UI
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLUI_h__
#define __OpenGLUI_h__

#include <memory>
#include "GeometryPrimitives.h"
#include "OpenGLObject.h"
#include "OpenGLViewer.h"
#include "OpenGLUiImgui.h"
#include "OpenGLWindow.h"

class OpenGLUI : public OpenGLObject
{public:typedef OpenGLObject Base;using Base::name;using Base::Update_Polygon_Mode;
	Box<2> box;
	int mouse_x=-1,mouse_y=-1;
	int left=0,right=0,middle=0,scroll=0;
	bool is_mouse_hovering_window=false;
	Vector2 pos=Vector2((real)0,(real)0);	////starting from the top-left corner
	OpenGLViewer* opengl_viewer=nullptr;

	OpenGLUI(){name="basic_ui";polygon_mode=PolygonMode::Fill;}

	virtual void Set_Pos(const Vector2& _pos){pos=_pos;}

	virtual void Display() const
    {
		Update_Polygon_Mode();
		ImguiGlut::New_Frame_Imgui(mouse_x,mouse_y,0,right,middle,scroll);
		ImGui::SetNextWindowPos(ImVec2((float)pos[0],(float)pos[1]));
		ImGui::SetNextWindowSize(ImVec2(256,256));
		ImGui::Begin("UI");
		ImGui::Text("Basic UI");
		ImGui::End();

		ImGui::Render();
    }

	virtual bool Mouse_Drag(int x,int y,int w,int h)
	{
		mouse_x=x;mouse_y=h-y;
		is_mouse_hovering_window=ImguiGlut::Update_Mouse_Position(mouse_x,mouse_y);
		return is_mouse_hovering_window;
	}

	virtual bool Mouse_Click(int _left,int _right,int _middle,int x,int y,int w,int h)
	{
		left=_left;
		right=_right;middle=_middle;mouse_x=x;mouse_y=h-y;
		is_mouse_hovering_window=ImguiGlut::Update_Mouse_Position(mouse_x,mouse_y);
		return is_mouse_hovering_window;
	}

	virtual bool Mouse_Scroll(int wheel,int direction,int x,int y)
	{scroll=direction;return is_mouse_hovering_window;}
};


static const int text_size=64;
static char input_text[text_size]={0};

class OpenGLUICommand : public OpenGLUI
{public:typedef OpenGLUI Base;using Base::name;Base::Update_Polygon_Mode;using Base::visible;
	Array<std::string> commands;
	OpenGLUICommand():Base(){name="ui_command";}

	bool toggle_exit=false;

	virtual void Display() const
    {
		if(!visible)return;

		Update_Polygon_Mode();
		ImguiGlut::New_Frame_Imgui(mouse_x,mouse_y,0,right,middle,scroll);

		float w=460;float h=42;
		ImGui::SetNextWindowSize(ImVec2(w,h));
		ImGui::SetNextWindowPos(ImVec2(Win_Width()*.5f-w*.5f,Win_Height()*0.9f));
		ImGui::Begin("Command Window",0,ImGuiWindowFlags_NoTitleBar|ImGuiWindowFlags_NoResize);
		ImGui::SetWindowFontScale(2.1f);
		ImGui::PushItemWidth(w-60);
		ImGui::PushStyleColor(ImGuiCol_Text,ImVec4(1.f,1.f,0.f,1.f));
		ImGui::PushStyleColor(ImGuiCol_PopupBg,ImVec4(.1f,1.f,.1f,1.f));
		ImGui::SetKeyboardFocusHere();
		bool is_done=ImGui::InputText("Cmd",input_text,text_size,ImGuiInputTextFlags_EnterReturnsTrue);
		ImGui::PopStyleColor(2);
		ImGui::PopItemWidth();
		ImGui::End();
		ImGui::Render();

		if(is_done){for(int i=0;i<text_size;i++){input_text[i]=0;}}
		if(toggle_exit){
			auto non_const_ptr=const_cast<OpenGLUICommand*>(this);
			non_const_ptr->visible=false;
			non_const_ptr->toggle_exit=false;}
    }

	////interactive command line
	virtual bool Keyboard(unsigned char key,int x,int y,bool down)
	{
		if(key==(unsigned char)13&&down&&!visible){visible=true;return false;}
		if(!visible){return false;}
		////Enter and get command
		if(key==(unsigned char)13&&down&&visible){
			std::string cmd=Get_Command();
			std::cout<<"Cmd: "<<cmd<<std::endl;
			commands.push_back(cmd);
			if(opengl_viewer!=nullptr)opengl_viewer->Toggle_Command(cmd);

			//if(cmd=="exit"){toggle_exit=true;}

			////default action: exit UI after a command input
			toggle_exit=true;
		}

		////Process input key
		ImguiGlut::Update_Keyboard(key,x,y,down);
		return true;
	}

	virtual bool Keyboard_Special(int key,int x,int y,bool down)
	{
		if(!visible){return false;}
		ImguiGlut::Update_Keyboard_Special(key,x,y,down);
		return true;
	}

	std::string Get_Command() const
	{
		std::string cmd;
		int i=0;for(;i<text_size;i++){if(input_text[i]==0)break;}
		cmd.append(&input_text[0],i);
		return cmd;
	}
};

class OpenGLUICurvePlot : public OpenGLUI
{public:typedef OpenGLUI Base;using Base::name;Base::Update_Polygon_Mode;using Base::visible;
	Array<Array<float> > values;
	Array<int> index;
	Array<std::string> line_names;
	Array<std::string> value_names;
	std::string plot_name="";

	OpenGLUICurvePlot(){name="ui_curve_plot";}

	virtual void Initialize(){}

	virtual void Display() const
	{
		if(!visible)return;

		Update_Polygon_Mode();
		ImguiGlut::New_Frame_Imgui(mouse_x,mouse_y,0,right,middle,scroll);
		ImGui::SetNextWindowPos(ImVec2((float)pos[0],(float)pos[1]));
		ImGui::SetNextWindowSize(ImVec2(196.f,(float)(64*(values.size()+1))));
		ImGui::Begin(plot_name.c_str());
		ImGui::Text(plot_name.c_str());
		for(int i=0;i<(int)values.size();i++){
			ImGui::PlotLines(line_names[i].c_str(),&values[i][0],(int)values.size(),index[i],value_names[i].c_str(),-1.f,1.f,ImVec2(0,64));}
		ImGui::End();

		ImGui::Render();
	}

	void Resize_Lines(int line_n,int value_n=256)
	{
		values.resize(line_n);
		index.resize(line_n);
		line_names.resize(line_n);
		value_names.resize(line_n);
		for(int i=0;i<line_n;i++){
			values[i].resize(value_n,0);}
	}

	void Next_Step()
	{
		for(int i=0;i<(int)values.size();i++){
			index[i]++;if(index[i]>=values.size())index[i]=0;}
	}
};
#endif
