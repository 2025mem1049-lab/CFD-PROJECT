## Developing Flow Between Two Parallel Plates
Developed a numerical CFD model for hydrodynamically developing flow between parallel plates using the ψ–ω formulation for Re = 100, 500 and 1000.
Investigated hydrodynamic entry length, boundary-layer growth, velocity-profile development, and Reynolds-number effects, comparing the results with the analytical Poiseuille solution.
Implemented upwind convection, central-difference diffusion, adaptive time stepping, and Gauss–Seidel Poisson solution, with residual-based convergence monitoring.
Demonstrated the dependence of entry length on Reynolds number and identified near-wall grid resolution as an important factor at higher Reynolds numbers.

## Flow Over a Flat Plate – Boundary Layer Analysis
Developed a 2-D incompressible CFD model for laminar flow over a flat plate using the Stream Function–Vorticity (ψ–ω) formulation.
Computed velocity and vorticity fields, streamlines, boundary-layer thickness, and velocity profiles, with results compared against the Blasius similarity solution.
Implemented numerical post-processing for boundary-layer thickness and velocity profiles at multiple streamwise locations.

## LID driven cavity
Developed a 2-D lid-driven cavity flow solver using the Stream Function–Vorticity formulation to analyze recirculating flow and vortex formation at Re = 100.
Implemented central finite differences, explicit Euler time integration, adaptive CFL-based time stepping, and Gauss–Seidel iteration for solving the governing equations.
Applied Thom’s wall-vorticity boundary condition and analyzed stream function, vorticity, velocity-field contours, and numerical convergence.
