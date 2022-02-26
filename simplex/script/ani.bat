set output=%1
set res=%2
set image=%3
set file=%4

if "%output%"=="" goto :error1
if "%res%"=="" set res=1024
rem if "%image%"=="" set image=_image
if "%file%"=="" set file=out

echo "Create simplex animation!"
echo %output%
echo %res%
echo %image%
echo %file%

ffmpeg -framerate 25 -i %output%/%image%/%%4d.png -c:v libx264 -vf "crop=trunc(iw/2)*2:trunc(ih/2)*2" -pix_fmt yuv420p %output%/%image%/%file%.mp4

goto :finish_ffmpeg_successfully

:error1
echo [Error]: output folder not specified!

:finish_ffmpeg_successfully