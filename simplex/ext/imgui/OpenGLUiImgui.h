// ImGui GLFW binding with OpenGL3 + shaders
// In this binding, ImTextureID is used to store an OpenGL 'GLuint' texture identifier. Read the FAQ about ImTextureID in imgui.cpp.

// You can copy and use unmodified imgui_impl_* files in your project. See main.cpp for an example of using this.
// If you use this binding you'll need to call 4 functions: ImGui_ImplXXXX_Init(), ImGui_ImplXXXX_NewFrame(), ImGui::Render() and ImGui_ImplXXXX_Shutdown().
// If you are new to ImGui, see examples/README.txt and documentation at the top of imgui.cpp.
// https://github.com/ocornut/imgui

//struct GLFWwindow;
//
//IMGUI_API bool        ImGui_ImplGlfwGL3_Init(GLFWwindow* window, bool install_callbacks);
//IMGUI_API void        ImGui_ImplGlfwGL3_Shutdown();
//IMGUI_API void        ImGui_ImplGlfwGL3_NewFrame();
//
//// Use if you want to reset your rendering device without losing ImGui state.
//IMGUI_API void        ImGui_ImplGlfwGL3_InvalidateDeviceObjects();
//IMGUI_API bool        ImGui_ImplGlfwGL3_CreateDeviceObjects();
//
//// GLFW callbacks (installed by default if you enable 'install_callbacks' during initialization)
//// Provided here if you want to chain callbacks.
//// You can also handle inputs yourself and use those as a reference.
//IMGUI_API void        ImGui_ImplGlfwGL3_MouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
//IMGUI_API void        ImGui_ImplGlfwGL3_ScrollCallback(GLFWwindow* window, double xoffset, double yoffset);
//IMGUI_API void        ImGui_ImplGlfwGL3_KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
//IMGUI_API void        ImGui_ImplGlfwGL3_CharCallback(GLFWwindow* window, unsigned int c);

#ifndef __OpenGLUiImgui_h__
#define __OpenGLUiImgui_h__
#include "imgui.h"

namespace ImguiGlut
{
static const int glut_special_key_map_offset=256;
IMGUI_API bool Initialize_Imgui();
IMGUI_API void New_Frame_Imgui(const int mouse_x,const int mouse_y,const int left,const int right,const int middle,const int scroll);
IMGUI_API void Shut_Down_Imgui();
IMGUI_API bool Update_Mouse_Position(const int mouse_x,const int mouse_y);
IMGUI_API bool Update_Keyboard(const unsigned char key,const int mouse_x,const int mouse_y,const bool down);
IMGUI_API bool Update_Keyboard_Special(const int key,const int mouse_x,const int mouse_y,const bool down);
//IMGUI_API void Clear_Text_Input();
};

#endif // !__OpenGLUiImgui_h__