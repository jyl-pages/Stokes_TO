# OpenGL Viewer

SimpleX implements an OpenGL viewer to help you visualize your simulation data. The viewer implementation relies on an old-fashioned OpenGL window, glut, but with the latest GLSL shading language support. The entire source code for the viewer lies in two folder: src/viewer and proj/opengl_viewer. You need to compile the executable for the viewer before using it. You may use an environment variable for a quick access (by default it is %o% specified in the install_env_var.bat script, but you can always customize your own).



## Guidelines for basic functions

The basic usage of the viewer is very straightforward. You can simply type the following command in a command line:
```
o -m mode_name -o folder_name
```

The mode name specifies the data mode (e.g, fluid, solid, topo, mesh, etc.). The folder name specifies where the data you want to visualize is stored (e.g., output). Then you can see an OpenGL window with the first frame displayed on the screen. The viewer pre-defined a set of hotkeys for the basic animation functions. These keys include:

* p: play/stop the animation
* [: previous frame
* ]: next frame
* r: get back to the first frame
* \> : scale up the length scale of the data (e.g., the length of velocity vectors, you need to refresh the frame by pressing ] or [ to see its update)
* < : scale down the length scale of the data
* w: turn on/off offscreen rendering
* q: quite the window
* K: print all the hotkeys in the command line



At the same time, each visualization mode contains a set of customized hotkeys. The key-function mapping of these hotkeys are usually defined in `proj/opengl_viewer/OpenGLViewerDriver.cpp` (typically as an input parameter for the function `Set_Visibility`). You may look them up and customize the key mappings by reading the source code file.

## Show Element Index in opengl_viewer

### For users of opengl_viewer

When debugging, you can use opengl_viewer to print the index of elements showing on screen.

This feature now support following OpenGLObject:

+ OpenGLTrackerCircles

Take OpenGLTrackerCircles as an example.

1. Execute opengl_viewer.exe with command-line, display some `OpenGLTrackerCircles` produced by your simulation.

<img width="50%" src="/_static/debug/debug_1.png"/>

2. Press Shift + O to enter debugging mode. The color of displayed elements will change.

<img width="50%" src="/_static/debug/debug_2.png"/>

3. Click to select the element you want to view. Its color will turn white.

<img width="50%" src="/_static/debug/debug_3.png"/>

4. The index number of the selected element will be printed in commane-lines (16 in this example).

<img width="50%" src="/_static/debug/debug_4.png"/>

5. Press Shift + O again to exit debug mode.

### For developers of opengl_viewer

This debug tool in opengl_viewer is implemented with color-coding technique. Which means to give every element a unique color depending on its index. On clicking, opengl_viewer will detect the color of the pixel you clicked, and decide the element of selection.

If you want to debug another object, which is a derived-class of `OpenGLObject` (say `SomeObject` for convenience), you need to do following things:

1. Inherit and implement `int Color_Code_Size()` in `SomeObject`. This function should return the number of elements in SomeObject. Note that it may be more than one, as in the case of `OpenGLTrackerCircles`, which is a set of circles.

1. When initializing `SomeObject`, add shader `vpos_model` with
    ```cpp
    Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vpos_model"));
    ```
    Which means no lightning, no ambient, no reflection... everything is displayed with just its color.

2. Modify `Display()` function in `SomeObject`. if its `shading_mode==ShadingMode::ColorEncoding`, display every element in `SomeObject` (it may contain more than one) with shader `vpos_model`. The color `color_code` of `i`-th element (indexing from 0) is given by
    ```cpp
    if (i == color_code_selection) color_code = OpenGLColor::White();
    else color_code = ID_To_Color(i + color_code_offset);
    ```
    Where `color_code_selection` and `color_code` are members of `OpenGLObject`. They are automatically handled by `OpenGLObject`, and `SomeObject` should not specifically deal with them.

1. (Optional) inherit and implement `Output_Debug_Info()` in `SomeObject`. Maybe you want to show other information, such as the number of elements, except for what is already printed in `OpenGLObject::Output_Debug_Info(out);`

You can refer to `class OpenGLTrackerCircles` for more details.


## Build your customized viewer mode

**TODO**



## Configurable viewer mode with json

**TODO**



## Generate animation

### Generate animation video with ffmpeg

The OpenGL viewer in SimpleX supports the automated generation of animations. For any simulation data you visualize in the viewer, press 'w' to enter the offscreen rendering mode. Then, after you start to play the animation, the viewer will write each frame into an image stored in the folder "[output]/_images." To get an animated video, you need to use ffmpeg to put the frames together. You need to have ffmpeg installed on your computer before recording animation.

1. In the viewer, press `w` to turn on the offscreen rendering mode. Then press `p` to play the animation in the viewer. A sequence of frames will be written to `[your output path]/_images` (e.g., `water/_images`)

2. Use ffmpeg to generate video.

    ffmpeg -framerate 25 -i [output]/_images/%4d.png -c:v libx264 -vf scale=640:-1 -pix_fmt yuv420p [output]/_images/out.mp4

(Reference: https://trac.ffmpeg.org/wiki/Slideshow)



### Use the %a% shortcut command

We have a pre-installed environment variable %a% (see the previous section for environment variable) to shorten the command of ffmpeg. The usage is very simple. Run the following command:

```
%a% output/_images
```

If everything runs successfully, you will get an .mp4 video composed of all the frames within the directory you specified in the command line.



### Generate animation video on Linux

If you use WSL, you may create a bash file for templaterized command call (see a bash file in `simplex/script/ani.bash` for example). You may also add command aliases in the WSL .bashrc file (like Windows environmental variables) for quick access:

    cd ~
    (goto home dir)
    nano .bashrc
    (open the .bashrc file with nano)

Then append the following text to the end of the file and save the file:

    # customized aliases
    alias ani='[Your path]/simplex/script/ani.bash'

    source .bashrc
    (update the aliases)



## Interactive visualization

**TODO**

# Authors

Bo Zhu, Mengdi Wang
