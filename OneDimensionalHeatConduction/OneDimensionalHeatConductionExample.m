% clear all
% close all
% clc

% Defing the necessary parameters
n = 500;
dx = 0.1*ones(n,1);
temp = 100*ones(n,1);
dt = 1;
tol = 1e-14;
max_iters = 1000;

% Defining the functions governing the necessary properties
function obj = k_f(obj)
    arguments
        obj (1,1) OneDimensionalHeatConductionSimulation
    end
    obj.k = 50*ones(obj.sim_size,1);
end

function obj = dcp_f(obj)
    arguments
        obj (1,1) OneDimensionalHeatConductionSimulation
    end
    obj.dcp = 50*ones(obj.sim_size,1);
end

function obj = q_f(obj)
    arguments
        obj (1,1) OneDimensionalHeatConductionSimulation
    end
    obj.q = 20*ones(obj.sim_size,1);
end

function obj = qt_f(obj)
    arguments
        obj (1,1) OneDimensionalHeatConductionSimulation
    end
    obj.qt = 0*ones(obj.sim_size,1);
end

% Creating the object
sim = OneDimensionalHeatConductionSimulation(n, dx, temp, dt, tol, max_iters, @k_f, @dcp_f, @q_f, @qt_f);

% Setting the boundary conditions
sim = sim.set_dirichlet_boundary_condition(0, 100);
sim = sim.set_robin_boundary_condition(1,200,250);

% Calculating the transient state and plotting the temperature
figure(1);
hold on
for i = 1:5
    sim = sim.unsteady(100);
    sim.plot_temperature();
end
hold off

% Plotting the residual history
figure(2);
sim.plot_residual_history();

sim = sim.delete_residual_history();

% Calculating the steady state and plotting the temperature
figure(3)
sim = sim.steady();
sim.plot_temperature()

figure(4);
sim.plot_residual_history();

