set proj_name=%1
set c=%2

cd %simplex%
git pull origin master
cd %complex%
git pull origin master

if not "%proj_name%"=="" (%b% %proj_name% %c%)
