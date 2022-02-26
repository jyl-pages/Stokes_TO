# Overview

<img width="100%" src="/_static/simplex_structure.jpg"/>

### Must read:

* **Eigen** is used as the vector math library. The most frequently used Eigen types include `Vector2, Vector3, Vector2i, Vector3i, VectorD and VectorDi` (within a dimensional templaterized class), MatrixD, etc.
* The C++ **std** library is used as the data containers, including `std::vector (Array), array (with fixed size, ArrayF), std::unordered_map (Hashtable), std::unordered_set (Hashset), std::unordered_multimap (HashtableMultiValue)`. The names within the parentheses are the aliased names of these containers used in SimpleX.
* All the Eigen, std::, and other data names being used frequently are defined collectively in `src/common/Common.h`. This header file is one of the most frequently included file in SimpleX.
* The codebase uses **double** floating-point precision (type **real**), though it can be converted to single precision by defining real as float in Common.h.  
* The codebase is written mostly in **C++ 11**, but it also covers 14 and 17 features occasionally (e.g., `constexpr` condition).  
* The codebase is **templaterized** as much as possible.
* The main application of templaterization is with its **dimension parameter d**, meaning that in most cases, for each simulator, you need to implement a single copy of code to support both 2D and 3D simulations.
*  The code uses the **right-hand** coordinate system.
* The Eigen **dense storage is column-major** ((e.g., MatrixD, MatrixX, VectorD) ).
* The Eigen **sparse storage is row-major** (E.g., SparseMatrixT). Attention that this convention is different from the Eigen default setting!
