Step 1
Run following command:
# python EULER-1A.py exam.pkl
This file EULER-1A.py describes a exam, running project simplex/fluid_euler, with 5 tests where resolution=50,100,150,200,250 correspondingly.

Step 2
Suppose your CPU have 4 cores, so you want 4 processes parallely running all tests (we can handle arbitrary number of proccesses and tests), so just parallely run 4 commands in shell at the same time:
# python test_runner.py exam.pkl 4 0
# python test_runner.py exam.pkl 4 1
# python test_runner.py exam.pkl 4 2
# python test_runner.py exam.pkl 4 3
Which means the 0,1,2,3 part of all 4 parts.

See EULER-1A.py to know how to construct a such file. You can create your own with a different name.