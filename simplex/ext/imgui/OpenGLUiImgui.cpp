// ImGui GLFW binding with OpenGL3 + shaders
// In this binding, ImTextureID is used to store an OpenGL 'GLuint' texture identifier. Read the FAQ about ImTextureID in imgui.cpp.

// You can copy and use unmodified imgui_impl_* files in your project. See main.cpp for an example of using this.
// If you use this binding you'll need to call 4 functions: ImGui_ImplXXXX_Init(), ImGui_ImplXXXX_NewFrame(), ImGui::Render() and ImGui_ImplXXXX_Shutdown().
// If you are new to ImGui, see examples/README.txt and documentation at the top of imgui.cpp.
// https://github.com/ocornut/imgui

#include <imgui.h>
#include <imgui_internal.h>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <iostream>
#include "OpenGLUiImgui.h"

namespace ImguiGlut
{
static double       g_Time = 0.0f;
static GLuint       g_FontTexture = 0;
static int          g_ShaderHandle = 0, g_VertHandle = 0, g_FragHandle = 0;
static int          g_AttribLocationTex = 0, g_AttribLocationProjMtx = 0;
static int          g_AttribLocationPosition = 0, g_AttribLocationUV = 0, g_AttribLocationColor = 0;
static unsigned int g_VboHandle = 0, g_VaoHandle = 0, g_ElementsHandle = 0;
static int          oldTimeSinceStart = 0;

//// This is the main rendering function that you have to implement and provide to ImGui (via setting up 'RenderDrawListsFn' in the ImGuiIO structure)
//// If text or lines are blurry when integrating ImGui in your engine:
//// - in your Render function, try translating your projection matrix by (0.5f,0.5f) or (0.375f,0.375f)
void ImGui_ImplGlut_RenderDrawLists(ImDrawData* draw_data)
{
    // Avoid rendering when minimized, scale coordinates for retina displays (screen coordinates != framebuffer coordinates)
    ImGuiIO& io = ImGui::GetIO();
    int fb_width = (int)(io.DisplaySize.x * io.DisplayFramebufferScale.x);
    int fb_height = (int)(io.DisplaySize.y * io.DisplayFramebufferScale.y);
    if (fb_width == 0 || fb_height == 0)
        return;
    draw_data->ScaleClipRects(io.DisplayFramebufferScale);

    // Backup GL state
    GLint last_active_texture; glGetIntegerv(GL_ACTIVE_TEXTURE, &last_active_texture);
    glActiveTexture(GL_TEXTURE0);
    GLint last_program; glGetIntegerv(GL_CURRENT_PROGRAM, &last_program);
    GLint last_texture; glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture);
    GLint last_array_buffer; glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &last_array_buffer);
    GLint last_element_array_buffer; glGetIntegerv(GL_ELEMENT_ARRAY_BUFFER_BINDING, &last_element_array_buffer);
    GLint last_vertex_array; glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &last_vertex_array);
    GLint last_blend_src_rgb; glGetIntegerv(GL_BLEND_SRC_RGB, &last_blend_src_rgb);
    GLint last_blend_dst_rgb; glGetIntegerv(GL_BLEND_DST_RGB, &last_blend_dst_rgb);
    GLint last_blend_src_alpha; glGetIntegerv(GL_BLEND_SRC_ALPHA, &last_blend_src_alpha);
    GLint last_blend_dst_alpha; glGetIntegerv(GL_BLEND_DST_ALPHA, &last_blend_dst_alpha);
    GLint last_blend_equation_rgb; glGetIntegerv(GL_BLEND_EQUATION_RGB, &last_blend_equation_rgb);
    GLint last_blend_equation_alpha; glGetIntegerv(GL_BLEND_EQUATION_ALPHA, &last_blend_equation_alpha);
    GLint last_viewport[4]; glGetIntegerv(GL_VIEWPORT, last_viewport);
    GLint last_scissor_box[4]; glGetIntegerv(GL_SCISSOR_BOX, last_scissor_box);
    GLboolean last_enable_blend = glIsEnabled(GL_BLEND);
    GLboolean last_enable_cull_face = glIsEnabled(GL_CULL_FACE);
    GLboolean last_enable_depth_test = glIsEnabled(GL_DEPTH_TEST);
    GLboolean last_enable_scissor_test = glIsEnabled(GL_SCISSOR_TEST);

    // Setup render state: alpha-blending enabled, no face culling, no depth testing, scissor enabled
    glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_SCISSOR_TEST);

    // Setup viewport, orthographic projection matrix
    glViewport(0, 0, (GLsizei)fb_width, (GLsizei)fb_height);
    const float ortho_projection[4][4] =
    {
        { 2.0f/io.DisplaySize.x, 0.0f,                   0.0f, 0.0f },
        { 0.0f,                  2.0f/-io.DisplaySize.y, 0.0f, 0.0f },
        { 0.0f,                  0.0f,                  -1.0f, 0.0f },
        {-1.0f,                  1.0f,                   0.0f, 1.0f },
    };
    glUseProgram(g_ShaderHandle);
    glUniform1i(g_AttribLocationTex, 0);
    glUniformMatrix4fv(g_AttribLocationProjMtx, 1, GL_FALSE, &ortho_projection[0][0]);
    glBindVertexArray(g_VaoHandle);

    for (int n = 0; n < draw_data->CmdListsCount; n++)
    {
        const ImDrawList* cmd_list = draw_data->CmdLists[n];
        const ImDrawIdx* idx_buffer_offset = 0;

        glBindBuffer(GL_ARRAY_BUFFER, g_VboHandle);
        glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)cmd_list->VtxBuffer.Size * sizeof(ImDrawVert), (const GLvoid*)cmd_list->VtxBuffer.Data, GL_STREAM_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g_ElementsHandle);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, (GLsizeiptr)cmd_list->IdxBuffer.Size * sizeof(ImDrawIdx), (const GLvoid*)cmd_list->IdxBuffer.Data, GL_STREAM_DRAW);

        for (int cmd_i = 0; cmd_i < cmd_list->CmdBuffer.Size; cmd_i++)
        {
            const ImDrawCmd* pcmd = &cmd_list->CmdBuffer[cmd_i];
            if (pcmd->UserCallback)
            {
                pcmd->UserCallback(cmd_list, pcmd);
            }
            else
            {
                glBindTexture(GL_TEXTURE_2D, (GLuint)(intptr_t)pcmd->TextureId);
                glScissor((int)pcmd->ClipRect.x, (int)(fb_height - pcmd->ClipRect.w), (int)(pcmd->ClipRect.z - pcmd->ClipRect.x), (int)(pcmd->ClipRect.w - pcmd->ClipRect.y));
                glDrawElements(GL_TRIANGLES, (GLsizei)pcmd->ElemCount, sizeof(ImDrawIdx) == 2 ? GL_UNSIGNED_SHORT : GL_UNSIGNED_INT, idx_buffer_offset);
            }
            idx_buffer_offset += pcmd->ElemCount;
        }
    }

    // Restore modified GL state
    glUseProgram(last_program);
    glBindTexture(GL_TEXTURE_2D, last_texture);
    glActiveTexture(last_active_texture);
    glBindVertexArray(last_vertex_array);
    glBindBuffer(GL_ARRAY_BUFFER, last_array_buffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, last_element_array_buffer);
    glBlendEquationSeparate(last_blend_equation_rgb, last_blend_equation_alpha);
    glBlendFuncSeparate(last_blend_src_rgb, last_blend_dst_rgb, last_blend_src_alpha, last_blend_dst_alpha);
    if (last_enable_blend) glEnable(GL_BLEND); else glDisable(GL_BLEND);
    if (last_enable_cull_face) glEnable(GL_CULL_FACE); else glDisable(GL_CULL_FACE);
    if (last_enable_depth_test) glEnable(GL_DEPTH_TEST); else glDisable(GL_DEPTH_TEST);
    if (last_enable_scissor_test) glEnable(GL_SCISSOR_TEST); else glDisable(GL_SCISSOR_TEST);
    glViewport(last_viewport[0], last_viewport[1], (GLsizei)last_viewport[2], (GLsizei)last_viewport[3]);
    glScissor(last_scissor_box[0], last_scissor_box[1], (GLsizei)last_scissor_box[2], (GLsizei)last_scissor_box[3]);
}

static const char* ImGui_ImplGlut_GetClipboardText(void* user_data)
{
    //return glfwGetClipboardString((GLFWwindow*)user_data);
	return nullptr;
}

static void ImGui_ImplGlut_SetClipboardText(void* user_data, const char* text)
{
    //glfwSetClipboardString((GLFWwindow*)user_data, text);
}

//void ImGui_ImplGlfwGL3_KeyCallback(GLFWwindow*, int key, int, int action, int mods)
//{
//    ImGuiIO& io = ImGui::GetIO();
//    if (action == GLFW_PRESS)
//        io.KeysDown[key] = true;
//    if (action == GLFW_RELEASE)
//        io.KeysDown[key] = false;
//
//    (void)mods; // Modifiers are not reliable across systems
//    io.KeyCtrl = io.KeysDown[GLFW_KEY_LEFT_CONTROL] || io.KeysDown[GLFW_KEY_RIGHT_CONTROL];
//    io.KeyShift = io.KeysDown[GLFW_KEY_LEFT_SHIFT] || io.KeysDown[GLFW_KEY_RIGHT_SHIFT];
//    io.KeyAlt = io.KeysDown[GLFW_KEY_LEFT_ALT] || io.KeysDown[GLFW_KEY_RIGHT_ALT];
//    io.KeySuper = io.KeysDown[GLFW_KEY_LEFT_SUPER] || io.KeysDown[GLFW_KEY_RIGHT_SUPER];
//}
//

void Key_Callback(int key)
{
    ImGuiIO& io = ImGui::GetIO();
    //if (action == GLFW_PRESS)
    //    io.KeysDown[key] = true;
    //if (action == GLFW_RELEASE)
    //    io.KeysDown[key] = false;

    //(void)mods; // Modifiers are not reliable across systems
    //io.KeyCtrl = io.KeysDown[GLFW_KEY_LEFT_CONTROL] || io.KeysDown[GLFW_KEY_RIGHT_CONTROL];
    //io.KeyShift = io.KeysDown[GLFW_KEY_LEFT_SHIFT] || io.KeysDown[GLFW_KEY_RIGHT_SHIFT];
    //io.KeyAlt = io.KeysDown[GLFW_KEY_LEFT_ALT] || io.KeysDown[GLFW_KEY_RIGHT_ALT];
    //io.KeySuper = io.KeysDown[GLFW_KEY_LEFT_SUPER] || io.KeysDown[GLFW_KEY_RIGHT_SUPER];
	std::cout<<"key "<<key<<std::endl;
}

void Char_Callback(unsigned int c)
{
    ImGuiIO& io = ImGui::GetIO();
    if (c > 0 && c < 0x10000)
        io.AddInputCharacter((unsigned short)c);
}

bool ImGui_ImplGlut_CreateFontsTexture()
{
    // Build texture atlas
    ImGuiIO& io = ImGui::GetIO();
    unsigned char* pixels;
    int width, height;
    io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);   // Load as RGBA 32-bits (75% of the memory is wasted, but default font is so small) because it is more likely to be compatible with user's existing shaders. If your ImTextureId represent a higher-level concept than just a GL texture id, consider calling GetTexDataAsAlpha8() instead to save on GPU memory.

    // Upload texture to graphics system
    GLint last_texture;
    glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture);
    glGenTextures(1, &g_FontTexture);
    glBindTexture(GL_TEXTURE_2D, g_FontTexture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

    // Store our identifier
    io.Fonts->TexID = (void *)(intptr_t)g_FontTexture;

    // Restore state
    glBindTexture(GL_TEXTURE_2D, last_texture);

    return true;
}

bool ImGui_ImplGlut_CreateDeviceObjects()
{
    // Backup GL state
    GLint last_texture, last_array_buffer, last_vertex_array;
    glGetIntegerv(GL_TEXTURE_BINDING_2D, &last_texture);
    glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &last_array_buffer);
    glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &last_vertex_array);

    const GLchar *vertex_shader =
        "#version 330\n"
        "uniform mat4 ProjMtx;\n"
        "in vec2 Position;\n"
        "in vec2 UV;\n"
        "in vec4 Color;\n"
        "out vec2 Frag_UV;\n"
        "out vec4 Frag_Color;\n"
        "void main()\n"
        "{\n"
        "	Frag_UV = UV;\n"
        "	Frag_Color = Color;\n"
        "	gl_Position = ProjMtx * vec4(Position.xy,0,1);\n"
        "}\n";

    const GLchar* fragment_shader =
        "#version 330\n"
        "uniform sampler2D Texture;\n"
        "in vec2 Frag_UV;\n"
        "in vec4 Frag_Color;\n"
        "out vec4 Out_Color;\n"
        "void main()\n"
        "{\n"
        "	Out_Color = Frag_Color * texture( Texture, Frag_UV.st);\n"
        "}\n";

    g_ShaderHandle = glCreateProgram();
    g_VertHandle = glCreateShader(GL_VERTEX_SHADER);
    g_FragHandle = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(g_VertHandle, 1, &vertex_shader, 0);
    glShaderSource(g_FragHandle, 1, &fragment_shader, 0);
    glCompileShader(g_VertHandle);
    glCompileShader(g_FragHandle);
    glAttachShader(g_ShaderHandle, g_VertHandle);
    glAttachShader(g_ShaderHandle, g_FragHandle);
    glLinkProgram(g_ShaderHandle);

    g_AttribLocationTex = glGetUniformLocation(g_ShaderHandle, "Texture");
    g_AttribLocationProjMtx = glGetUniformLocation(g_ShaderHandle, "ProjMtx");
    g_AttribLocationPosition = glGetAttribLocation(g_ShaderHandle, "Position");
    g_AttribLocationUV = glGetAttribLocation(g_ShaderHandle, "UV");
    g_AttribLocationColor = glGetAttribLocation(g_ShaderHandle, "Color");

    glGenBuffers(1, &g_VboHandle);
    glGenBuffers(1, &g_ElementsHandle);

    glGenVertexArrays(1, &g_VaoHandle);
    glBindVertexArray(g_VaoHandle);
    glBindBuffer(GL_ARRAY_BUFFER, g_VboHandle);
    glEnableVertexAttribArray(g_AttribLocationPosition);
    glEnableVertexAttribArray(g_AttribLocationUV);
    glEnableVertexAttribArray(g_AttribLocationColor);

#define OFFSETOF(TYPE, ELEMENT) ((size_t)&(((TYPE *)0)->ELEMENT))
    glVertexAttribPointer(g_AttribLocationPosition, 2, GL_FLOAT, GL_FALSE, sizeof(ImDrawVert), (GLvoid*)OFFSETOF(ImDrawVert, pos));
    glVertexAttribPointer(g_AttribLocationUV, 2, GL_FLOAT, GL_FALSE, sizeof(ImDrawVert), (GLvoid*)OFFSETOF(ImDrawVert, uv));
    glVertexAttribPointer(g_AttribLocationColor, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(ImDrawVert), (GLvoid*)OFFSETOF(ImDrawVert, col));
#undef OFFSETOF

    ImGui_ImplGlut_CreateFontsTexture();

    // Restore modified GL state
    glBindTexture(GL_TEXTURE_2D, last_texture);
    glBindBuffer(GL_ARRAY_BUFFER, last_array_buffer);
    glBindVertexArray(last_vertex_array);

    return true;
}

void ImGui_ImplGlut_InvalidateDeviceObjects()
{
    if (g_VaoHandle) glDeleteVertexArrays(1, &g_VaoHandle);
    if (g_VboHandle) glDeleteBuffers(1, &g_VboHandle);
    if (g_ElementsHandle) glDeleteBuffers(1, &g_ElementsHandle);
    g_VaoHandle = g_VboHandle = g_ElementsHandle = 0;

    if (g_ShaderHandle && g_VertHandle) glDetachShader(g_ShaderHandle, g_VertHandle);
    if (g_VertHandle) glDeleteShader(g_VertHandle);
    g_VertHandle = 0;

    if (g_ShaderHandle && g_FragHandle) glDetachShader(g_ShaderHandle, g_FragHandle);
    if (g_FragHandle) glDeleteShader(g_FragHandle);
    g_FragHandle = 0;

    if (g_ShaderHandle) glDeleteProgram(g_ShaderHandle);
    g_ShaderHandle = 0;

    if (g_FontTexture)
    {
        glDeleteTextures(1, &g_FontTexture);
        ImGui::GetIO().Fonts->TexID = 0;
        g_FontTexture = 0;
    }
}

//bool ImGui_ImplGlfwGL3_Init(GLFWwindow* window, bool install_callbacks)
bool Initialize_Imgui()
{
//    //g_Window = window;
//
//    ImGuiIO& io = ImGui::GetIO();
//    io.KeyMap[ImGuiKey_Tab] = GLFW_KEY_TAB;                         // Keyboard mapping. ImGui will use those indices to peek into the io.KeyDown[] array.
//    io.KeyMap[ImGuiKey_LeftArrow] = GLFW_KEY_LEFT;
//    io.KeyMap[ImGuiKey_RightArrow] = GLFW_KEY_RIGHT;
//    io.KeyMap[ImGuiKey_UpArrow] = GLFW_KEY_UP;
//    io.KeyMap[ImGuiKey_DownArrow] = GLFW_KEY_DOWN;
//    io.KeyMap[ImGuiKey_PageUp] = GLFW_KEY_PAGE_UP;
//    io.KeyMap[ImGuiKey_PageDown] = GLFW_KEY_PAGE_DOWN;
//    io.KeyMap[ImGuiKey_Home] = GLFW_KEY_HOME;
//    io.KeyMap[ImGuiKey_End] = GLFW_KEY_END;
//    io.KeyMap[ImGuiKey_Delete] = GLFW_KEY_DELETE;
//    io.KeyMap[ImGuiKey_Backspace] = GLFW_KEY_BACKSPACE;
//    io.KeyMap[ImGuiKey_Enter] = GLFW_KEY_ENTER;
//    io.KeyMap[ImGuiKey_Escape] = GLFW_KEY_ESCAPE;
//    io.KeyMap[ImGuiKey_A] = GLFW_KEY_A;
//    io.KeyMap[ImGuiKey_C] = GLFW_KEY_C;
//    io.KeyMap[ImGuiKey_V] = GLFW_KEY_V;
//    io.KeyMap[ImGuiKey_X] = GLFW_KEY_X;
//    io.KeyMap[ImGuiKey_Y] = GLFW_KEY_Y;
//    io.KeyMap[ImGuiKey_Z] = GLFW_KEY_Z;
//
//    io.RenderDrawListsFn = ImGui_ImplGlfwGL3_RenderDrawLists;       // Alternatively you can set this to NULL and call ImGui::GetDrawData() after ImGui::Render() to get the same ImDrawData pointer.
//    io.SetClipboardTextFn = ImGui_ImplGlfwGL3_SetClipboardText;
//    io.GetClipboardTextFn = ImGui_ImplGlfwGL3_GetClipboardText;
//    io.ClipboardUserData = g_Window;
//#ifdef _WIN32
//    io.ImeWindowHandle = glfwGetWin32Window(g_Window);
//#endif
//
//    //if (install_callbacks)
//    //{
//    //    glfwSetMouseButtonCallback(window, ImGui_ImplGlfwGL3_MouseButtonCallback);
//    //    glfwSetScrollCallback(window, ImGui_ImplGlfwGL3_ScrollCallback);
//    //    glfwSetKeyCallback(window, ImGui_ImplGlfwGL3_KeyCallback);
//    //    glfwSetCharCallback(window, ImGui_ImplGlfwGL3_CharCallback);
//    //}
//
//    return true;

    ImGuiIO& io = ImGui::GetIO();

    io.RenderDrawListsFn = ImGui_ImplGlut_RenderDrawLists;   // Alternatively you can set this to NULL and call ImGui::GetDrawData() after ImGui::Render() to get the same ImDrawData pointer.
    io.SetClipboardTextFn = ImGui_ImplGlut_SetClipboardText;
    io.GetClipboardTextFn = ImGui_ImplGlut_GetClipboardText;

    io.KeyMap[ImGuiKey_Tab]        = 9;                     // tab
    io.KeyMap[ImGuiKey_LeftArrow]  = GLUT_KEY_LEFT+glut_special_key_map_offset;         // Left
    io.KeyMap[ImGuiKey_RightArrow] = GLUT_KEY_RIGHT+glut_special_key_map_offset;        // Right
    io.KeyMap[ImGuiKey_UpArrow]    = GLUT_KEY_UP+glut_special_key_map_offset;           // Up
    io.KeyMap[ImGuiKey_DownArrow]  = GLUT_KEY_DOWN+glut_special_key_map_offset;         // Down
    io.KeyMap[ImGuiKey_PageUp]     = GLUT_KEY_PAGE_UP+glut_special_key_map_offset;      // PageUp
    io.KeyMap[ImGuiKey_PageDown]   = GLUT_KEY_PAGE_DOWN+glut_special_key_map_offset;    // PageDown
    io.KeyMap[ImGuiKey_Home]       = GLUT_KEY_HOME+glut_special_key_map_offset;         // Home
    io.KeyMap[ImGuiKey_End]        = GLUT_KEY_END+glut_special_key_map_offset;          // End
    io.KeyMap[ImGuiKey_Delete]     = 127;                   // Delete
    io.KeyMap[ImGuiKey_Backspace]  = 8;                     // Backspace
    io.KeyMap[ImGuiKey_Enter]      = 13;                    // Enter
    io.KeyMap[ImGuiKey_Escape]     = 27;                    // Escape
    io.KeyMap[ImGuiKey_A]          = 1;                     // ctrl-A
    io.KeyMap[ImGuiKey_C]          = 3;                     // ctrl-C
    io.KeyMap[ImGuiKey_V]          = 22;                    // ctrl-V
    io.KeyMap[ImGuiKey_X]          = 24;                    // ctrl-X
    io.KeyMap[ImGuiKey_Y]          = 25;                    // ctrl-Y
    io.KeyMap[ImGuiKey_Z]          = 26;                    // ctrl-Z

    return true;
}

void Shut_Down_Imgui()
{
    ImGui_ImplGlut_InvalidateDeviceObjects();
    ImGui::Shutdown();
}

IMGUI_API bool Update_Mouse_Position(const int mouse_x, const int mouse_y)
{
	ImGuiIO& io = ImGui::GetIO();io.MousePos = ImVec2((float)mouse_x, (float)mouse_y);
	return ImGui::IsPosHoveringAnyWindow(io.MousePos);
}

IMGUI_API bool Update_Keyboard(const unsigned char key,const int mouse_x,const int mouse_y,const bool down)
{
	ImGuiIO& io = ImGui::GetIO();
	if(down){
		//std::cout<<"Key = "<<key<<", "<<(unsigned short)key<<std::endl;
		io.KeysDown[(int)key]=true;
		io.AddInputCharacter(key);
	}
	else io.KeysDown[(int)key]=false;

	return false;
}


IMGUI_API bool Update_Keyboard_Special(const int key,const int mouse_x,const int mouse_y,const bool down)
{
	ImGuiIO& io = ImGui::GetIO();
	if(down){
		//std::cout<<"Spe Key down = "<<key<<", "<<key<<std::endl;
		io.KeysDown[key+glut_special_key_map_offset]=true;
	}
	else{
		//std::cout<<"Spe Key up = "<<key<<", "<<key<<std::endl;
		io.KeysDown[key+glut_special_key_map_offset]=false;
	}

	return false;
}

//IMGUI_API void Clear_Text_Input()
//{
//	std::cout<<"Clear text input"<<std::endl;
//	ImGuiContext& g=*ImGui::GetCurrentContext();
//	std::cout<<"Size: "<<g.InputTextState.Text.size()<<std::endl;
//	for(int i=0;i<g.InputTextState.Text.size();i++){
//		g.InputTextState.TempTextBuffer[i]='\0';
//	    g.InputTextState.Text[i]='\0';
//		g.InputTextState.InitialText[i]='\0';
//	}
//	g.InputTextState.ClearSelection();
//	ImGuiIO& io = ImGui::GetIO();
//	io.ClearInputCharacters();
//}

void New_Frame_Imgui(const int mouse_x, const int mouse_y, const int left, const int right, const int middle, const int scroll)
{
    if (!g_FontTexture)
        ImGui_ImplGlut_CreateDeviceObjects();

    ImGuiIO& io = ImGui::GetIO();

    //// Setup display size (every frame to accommodate for window resizing)
    int w, h;
    int display_w, display_h;

    //glfwGetWindowSize(g_Window, &w, &h);
    //glfwGetFramebufferSize(g_Window, &display_w, &display_h);
	w=glutGet(GLUT_WINDOW_WIDTH);h=glutGet(GLUT_WINDOW_HEIGHT);
	display_w=w;display_h=h;
    io.DisplaySize = ImVec2((float)w, (float)h);
    io.DisplayFramebufferScale = ImVec2(w > 0 ? ((float)display_w / w) : 0, h > 0 ? ((float)display_h / h) : 0);

    //// Setup time step
    //double current_time =  glfwGetTime();
    //io.DeltaTime = g_Time > 0.0 ? (float)(current_time - g_Time) : (float)(1.0f/60.0f);
    //g_Time = current_time;
	//// Setup time step
    int timeSinceStart = glutGet(GLUT_ELAPSED_TIME);
    float deltaTime = (timeSinceStart - oldTimeSinceStart)/1000.0f;
    io.DeltaTime = deltaTime > 0.0 ? deltaTime : (float)(1.0f / 60.0f);
    oldTimeSinceStart = timeSinceStart;

	//std::cout<<"mouse: "<<mouse_x<<", "<<mouse_y<<", "<<left<<", "<<right<<", "<<middle<<std::endl;

    //// Setup inputs
    // (we already got mouse wheel, keyboard keys & characters from glfw callbacks polled in glfwPollEvents())
    //if (glfwGetWindowAttrib(g_Window, GLFW_FOCUSED))
    {
        //glfwGetCursorPos(g_Window, &mouse_x, &mouse_y);
        io.MousePos = ImVec2((float)mouse_x, (float)mouse_y);   // Mouse position in screen coordinates (set to -1,-1 if no mouse / on another screen, etc.)
    }
    //else
    //{
    //    io.MousePos = ImVec2(-1,-1);
    //}

	io.MouseDown[0]=left>0;
	io.MouseDown[1]=right>0;
	io.MouseDown[2]=middle>0;

    //for (int i = 0; i < 3; i++)
    //{
    //    io.MouseDown[i] = g_MousePressed[i] || glfwGetMouseButton(g_Window, i) != 0;    // If a mouse press event came, always pass it as "mouse held this frame", so we don't miss click-release events that are shorter than 1 frame.
    //    g_MousePressed[i] = false;
    //}

    //io.MouseWheel = g_MouseWheel;
    //g_MouseWheel = 0.0f;
	io.MouseWheel = (float)scroll;
	//std::cout<<"mouse_wheel: "<<scroll<<", "<<io.MouseWheel<<std::endl;

    //// Hide OS mouse cursor if ImGui is drawing it
    //glfwSetInputMode(g_Window, GLFW_CURSOR, io.MouseDrawCursor ? GLFW_CURSOR_HIDDEN : GLFW_CURSOR_NORMAL);

    // Start the frame
    ImGui::NewFrame();
}
};
