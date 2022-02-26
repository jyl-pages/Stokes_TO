# Physics

Here is a very brief inventory list of what SimpleX can do currently regarding physics simulation. All the source files for physics simulation are located within `src/physics`. In general, each file implements an independent simulator (e.g., Euler fluid) or functional module (e.g., advection, projection). 



### Euler fluid:

The grid-based simulator is a standard one as elaborated in Bridson's book, composed of the MAC grid structure, the semi-Lagrangian advection, and the incompressibility projection. Its variations include the level-set free-surface simulator, surface tension simulator, alpha values for non-axis-aligned boundaries, and PIC/FLIP. Viscosity was also implemented as an independent file. The essential part for an Euler fluid simulator lies in its Poisson solver, facilitated with a parallel multigrid solver on grid, to accelerate the computation.

The GPU implementation of the Euler solver is on its way.   



### Particle fluid:

SPH and MPM were implemented for particle fluid simulation. The sample projects are in `proj/fluid_sph` and `proj/mpm`. 



### Soft body:

SimpleX implements the mass-spring model, the finite-element model, and position-based dynamics for soft body simulations. The sample projects include `proj/mass_spring`, `proj/soft_body_fem`, and `proj/soft_body_position`. 



### Rigid body:

2D and 3D rigid body dynamics were implemented in `physics/Rigid.h`. The project is`proj/rigid_body`.



### Coupling:

Immersed boundary method. 