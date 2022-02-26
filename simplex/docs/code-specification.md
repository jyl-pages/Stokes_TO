# Code specification and naming convention

### Naming convention (must read!):

* Class name: Combination of **first-letter-capitalized words, no underscore in between**, e.g., 
  ```c++
  class FluidEuler{...};
  ```
  
* Variable name: Combination of **all-lower-case words, connected by underscores**, e.g.,
  ```c++
  FluidEuler grid_fluid_example;
  ```
  
* Function name: Combination of **first-letter-capitalized words, connected by underscores**, e.g.,

* Namespace name: the same as class name.

* Headings: put the code author name in the heading of each file. If you are modifying an existing copy, add your name after others. E.g., 
  ```c++
//////////////////////////////////////////////////////////////////////////
// SPH fluid
// Copyright (c) (2018-), Xiangxin Kong, Yaorui Zhang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
  ```


* Header macro: `__ClassName_h__`, e.g.,

  ```c++
  #ifndef __FaceField_h__
  #define __FaceField_h__
  ...
  #endif
  ```


### Other important notes:

* Spacing: I typically use extremely tight spacing for code, by removing all the spaces and line breaks as much as I can, which might bother many first-time readers. This was a convention from the PhysBAM library. I won't ask you to follow this convention in your code. But we have a python script to re-format the mature code once in a while to unify the spacing style (by automatically removing the spaces).
* **Important:** If you are a Microsoft Visual Studio user, please turn off its auto format option (**TODO: add how**), to avoid reformatting the existing code in the repository. 
* Do not check in debug code (e.g., `std::cout<<"why this code doesn't work?!"<<std::endl`) into the library in the `\scr\` folder. 
* Be very careful when making significant changes to the existing library code (the same, meaning any code that is outside your own project folder). Ideally, you'd better let the original author know about your changes, because this might affect many other projects that you were not aware of.
* Add comments to the code as much as you can. E.g., if you re-implement some equations from a paper, put the name of the paper in your comment; if you re-implement some algorithms from StackOverFlow, add the link in your comment.  