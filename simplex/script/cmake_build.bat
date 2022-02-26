set proj_name=%1
set c=%2

if "%c%"=="c" (set root_path=%complex%) else (set root_path=%simplex%)
if "%c%"=="r" (set root_path=%reflex%)

set src_path=%root_path%\proj\%proj_name%
set build_path=%root_path%\build\%proj_name%\cmake-build-win
set bin_path=%root_path%\bin\win\%proj_name%\Release
if not exist %build_path% mkdir %build_path%

cd %build_path%

if exist CMakeCache.txt rm CMakeCache.txt

cmake --clean-first . -G "Visual Studio 15 2017 Win64" --config Release --target %src_path%
if not %ERRORLEVEL% == 0 goto :error1

start /b %proj_name%.sln

if not exist %bin_path% (mkdir %bin_path%)
cd %bin_path%

goto :finish_cmake_build_successfully

:error1
echo [Error]: cmake failed!

:finish_cmake_build_successfully