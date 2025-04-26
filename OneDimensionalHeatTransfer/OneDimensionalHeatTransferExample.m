% clear all
% close all
% clc

% Defing the necessary parameters
% Parameters are constant to allow for comparison between the simulation
% results and the analytical results
n = 51;                     % Number of control volumes
L = 400;                    % Length of the domain
dx = (L/(n-1))*ones(n,1);   % Length of each of the control volumes

% The number of nodes in the denominator is reduced by one due to the
% half-control-volume formulation of the boundary conditions used in
% "Patankar, S. Numerical Heat Transfer and Fluid Flow"

% This formulation will be changed later so that boundary conditions will
% be given for the boundary faces

dt = 0.5;

cp = 4200;                  % Heat capacity
k = 20;                     % Thermal conductivity
rho = 1000;                 % Density
q = 0.005;                  % Constant heat generation term
v = 0.0000001;              % Velocity
TR = 50;                    % Temperature of the right boundary

qt = 0;                     % Linear heat generation term


temp = 100*ones(n,1);       % Initial temperature field

tol = 1e-14;                % Residue tolerance
max_iters = 1000;

% Defining the functions governing the necessary properties
function obj = k_f(obj, k)
    obj.k = k*ones(obj.sim_size,1);
end

function obj = rho_f(obj, rho)
    obj.rho = rho*ones(obj.sim_size,1);
end

function obj = q_f(obj,q)
    obj.q = q*ones(obj.sim_size,1);
end

function obj = qt_f(obj, qt)
    obj.qt = qt*ones(obj.sim_size,1);
end

function obj = v_f(obj, v)
    obj.v = v*ones(obj.sim_size,1);
end

% Creating the object
sim = OneDimensionalHeatTransferSimulation(n, dx, temp, cp, dt, tol, max_iters, ...
    @(obj) k_f(obj,k), @(obj) rho_f(obj,rho), @(obj) q_f(obj, q), @(obj) qt_f(obj, qt), @(obj) v_f(obj, v));

% Setting the boundary conditions
sim = sim.set_dirichlet_boundary_condition(0, 0);
sim = sim.set_dirichlet_boundary_condition(1, TR);

% Calculating the steady state and plotting the temperature
figure(1);
sim = sim.steady();
sim.plot_temperature()

% Plotting the residual history
figure(2);
sim.plot_residual_history();


% Analytical solution for the following boundary conditions:
% Constant temperatures at the boundaries, temperature of the left boundary
% equal to one
function x = analytical_solution(x, k, rho, cp, q, L, TR, v)
    b=q/(rho*cp);
    a=k/(rho*cp);
    x=(b*(-L*exp(v*x/a)+x*(exp(L*v/a)-1)+L)+TR*v*(exp(v*x/a)-1))/(v*(exp(L*v/a)-1));
end

figure(1);
hold on
plot(0:1:L,arrayfun(@(x) analytical_solution(x, k, rho, cp, q, L, TR, v),0:1:L),"Color","r");
hold off

