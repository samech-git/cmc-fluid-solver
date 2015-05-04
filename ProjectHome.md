This project has started in Moscow State University in the department of Computational Mathematics and Cybernetics, math-physics sub-department.

Our team:

| **name** | **status** | **roles** |
|:---------|:-----------|:----------|
| Nikolai Sakharnykh | post-grad | lead, development |
| Nikolay Markovskiy | PhD | multi-GPU development |
| Sergey Berezin | PhD | scientific advisor |
| Vilen M. Paskonov | Dr. of Sci. | scientific advisor |

We're working on a program product for simulation of incompressible viscid fluid in 2D/3D domains with dynamic boundaries. We use full system of Navier-Stokes equations that includes continuity, impulse and energy equations. For numerical methods we use finite difference methods on regular grids. The plans are to support some explicit and half-implicit (based on ADI scheme) finite diffrence methods.

The GPU acceleration will be used to speed-up the solver in large 2D/3D domains. Targeting on NVIDIA hardware, CUDA platform.

Target applications are aerodynamic flows, blood simulation in heart and vessels and underwater flows simulation in seas.

More information later.