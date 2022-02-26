#!/bin/bash
$ chmod +x ani.bash
echo "Create simplex animation!"
echo "Param 1-folder: $1"
echo "Param 2-resolution: $1"
ffmpeg -framerate 50 -i $1/_images/%4d.png -c:v libx264 -vf scale=$2:-1 -pix_fmt yuv420p $1/_images/_video.mp4