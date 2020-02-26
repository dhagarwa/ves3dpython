% MAIN CLASSES  /  Functions
calc.m % main class, abstract, calc2, and calc3 derived from it. 
calc2.m  % 2D class for FFT function approximations
calc3.m  % 3D class for FFT function approximations

rk2.m  % naive RK2 ode solver (same interface as matlab's ode45) but without time adaptivity
rk4.m  % naive RK4 ode solver (same interface as matlab's ode45) but without time adaptivity

% EXAMPLES
poisson.m % any D poisson solver
stokes2.m % 2D stokes solver example.

advection3.m % example of advection-diffusion solver in 3D
advection2.m % example of 2D advection-diffusion solver

% TESTS
test_grad.m % tests gradient calculations

% Incompressible fluids
euler.m % Euler solver
monitor.m  % used in euler (incompressible, inviscid fluid)
calc3_functions.m % some functions for euler solver

euler2.m % 2D euler equations solver (incompressible, inviscid fluid)
monitor2.m % function used in euler2 solver







