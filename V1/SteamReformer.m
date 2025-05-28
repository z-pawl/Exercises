%% Parametry symulacji
length = 0.3;
radius = 0.05;

porosity = 0.7;
permeability = 1e-7;
inertia_coefficient = 0.088;

solid_phase_thermal_conductivity = 20.0;
catalyst_density = 2.5e5;

temp_inlet = 800;
v_inlet = 0.001;
temp_wall = 800;
p_atm = 101325;
SC = 2;

a = 0.89;               % Order of reaction with respect to methane
b = 0.05;               % Order of reaction with respect to water
A_st = 1.354e-3;        % Arrhenius constant
E_a = 122500;           % Activation energy
delta_G = -28.6e3;      % Change of standard Gibbs free energy of water-gas-shift reaction
delta_H_st = 206e3;     % Enthalpy change accompanied with methane/steam reforming reaction
delta_H_sh = -41.2e3;   % Enthalpy change accompanied with water-gas-shift reaction

% Stała, której użycie wynika z użycia we wzorach na parametry jako
% argumentu temperatury podzielonej przez tysiąc
K = 1/1000;
% Stała w celu konwersji jednostek z uP na Pa*s
uP_to_Pas = 10^(-7);

% Parametry gazów
% Wodór
M_H2 = 2.016;
therm_cond_H2 = 0.01 * [-1.8954*K^6 11.972*K^5 -31.939*K^4 47.763*K^3 -47.190*K^2 62.892*K 1.5040];
visc_H2 = uP_to_Pas*[-9.9892*K^6 62.966*K^5 -167.51*K^4 249.41*K^3 -244.34*K^2 299.78*K 15.553];
heat_cap_H2 = [-6.4725*K^6 46.903*K^5 -136.15*K^4 199.29*K^3 -150.55*K^2 56.036*K 21.157];
diffusion_volume_H2 = 6.12;

% Tlenek węgla
M_CO = 28.010;
therm_cond_CO = 0.01 * [-2.3224*K^6 13.379*K^5 -30.818*K^4 36.018*K^3 -23.186*K^2 13.999*K -0.2815];
visc_CO = uP_to_Pas*[-32.298*K^6 208.42*K^5 -572.14*K^4 883.75*K^3 875.90*K^2 793.65*K -4.9137];
heat_cap_CO = [-7.6538*K^6 37.756*K^5 -66.346*K^4 41.974*K^3 5.2062*K^2 -8.1781*K 30.429];
diffusion_volume_CO = 18.0;

% Dwutlenek węgla
M_CO2 = 44.009;
therm_cond_CO2 = 0.01 * [18.698*K^6 -101.12*K^5 216.83*K^4 -233.29*K^3 129.65*K^2 -27.018*K 2.8888];
visc_CO2 = uP_to_Pas*[-0.4564*K^6 14.450*K^5 -85.929*K^4 244.22*K^3 -432.49*K^2 680.07*K -20.434];
heat_cap_CO2 = [-35.992*K^6 214.58*K^5 -519.9*K^4 657.88*K^3 -471.33*K^2 204.60*K 4.3669];
diffusion_volume_CO2 = 26.7;

% Metan
M_CH4 = 16.043;
therm_cond_CH4 = 0.01 * [3.2774*K^6 -17.283*K^5 38.251*K^4 -47.440*K^3 37.413*K^2 1.8732*K 0.4796];
visc_CH4 = uP_to_Pas*[-22.920*K^6 140.48*K^5 -367.06*K^4 548.11*K^3 -543.82*K^2 529.37*K -9.9989];
heat_cap_CH4 = [61.321*K^6 -358.75*K^5 856.93*K^4 -1068.7*K^3 712.55*K^2 -178.59*K 47.964];
diffusion_volume_CH4 = 25.14;

% Para wodna
M_H2O = 18.015;
therm_cond_H2O = 0.01 * [4.1531*K^6 -18.974*K^5 35.993*K^4 -41.390*K^3 35.922*K^2 -7.9139*K 2.0103];
visc_H2O = uP_to_Pas*[19.591*K^6 -126.96*K^5 348.12*K^4 -522.38*K^3 419.50*K^2 244.93*K -6.7541];
heat_cap_H2O = [14.015*K^6 -79.409*K^5 181.54*K^4 -217.08*K^3 146.01*K^2 -41.205*K 37.373];
diffusion_volume_H2O = 13.1;

% Parametry solvera

% Tolerancje
tol_T = 1e-7;
tol_v = 1e-7;
tol_p = 1e-7;
tol_species = 1e-7;

% Współczynniki relaksacji
rel_T = 1;
rel_v = 0.7;
rel_p = 0.25;
rel_species = 0.2;

% Maksymalna liczba iteracji
max_iters = 10000;

%% Definicja siatki
sz = [18 3];

dx = length*ones(sz(1),1)/sz(1);
dr = radius*ones(1,sz(2))/sz(2);

% Offset cannot be equal to zero as that would cause division by zero errors
offset = 1e-15;
grid = Grid2D(dx,dr,offset);

%% Definicja komponentów
sim = Simulation2D;

vel_component = Velocity2DSIMPLE(sim, grid, sz, zeros(sz+[1 0]), zeros(sz+[0 1]), p_atm*ones(sz), porosity, @(x) ones(sz), @(x) ones(sz), @(x) 0, @(x) 0, @(x) 0, @(x) 0, tol_v, tol_p, rel_v, rel_p);
temp_component = TemperatureCondConv2DSteady(sim, grid, vel_component, sz, (temp_wall+temp_inlet)/2*ones(sz), @(x) ones(sz), @(x) ones(sz), @(x) zeros(sz), @(x) zeros(sz), tol_T, rel_T);

hydrogen = Species(sim, grid, vel_component, sz, zeros(sz), heat_cap_H2, visc_H2, therm_cond_H2, M_H2, diffusion_volume_H2, tol_species, rel_species);
carbon_monoxide = Species(sim, grid, vel_component, sz, zeros(sz), heat_cap_CO, visc_CO, therm_cond_CO, M_CO, diffusion_volume_CO, tol_species, rel_species);
carbon_dioxide = Species(sim, grid, vel_component, sz, zeros(sz), heat_cap_CO2, visc_CO2, therm_cond_CO2, M_CO2, diffusion_volume_CO2, tol_species, rel_species);
methane = Species(sim, grid, vel_component, sz, M_CH4/(M_CH4+SC*M_H2O)*ones(sz), heat_cap_CH4, visc_CO, therm_cond_CH4, M_CH4, diffusion_volume_CH4, tol_species, rel_species);
water = Species(sim, grid, vel_component, sz, SC*M_H2O/(M_CH4+SC*M_H2O)*ones(sz), heat_cap_H2O, visc_H2O, therm_cond_H2O, M_H2O, diffusion_volume_H2O, tol_species, rel_species);

reactions = ChemicalReactions(sim, grid, vel_component, temp_component, sz,a,b,A_st,E_a,delta_G,SC,catalyst_density,v_inlet*pi*radius^2/(1+SC),length*pi*radius^2, ...
    ["hydrogen" "carbon monoxide" "carbon dioxide" "methane" "water"], ...
    [hydrogen carbon_monoxide carbon_dioxide methane water]);

%% Definicja członów źródłowych
% Źródła komponentów prędkości
function src_lin_vx = linear_source_vx(vel, permeability, inertia_coefficient)
    % Interpolating the radial velocity components, the density and the 
    % viscosity to the position of the axial velocity components
    vr_int = Utilities.interpolate_center2faces((vel.vr(:,1:end-1)+vel.vr(:,2:end))/2,vel.grid.dx,1);
    rho_int = Utilities.interpolate_center2faces(vel.rho,vel.grid.dx,1);
    visc_int = Utilities.interpolate_center2faces(vel.visc,vel.grid.dx,1);

    % Source term = vx*(-μ/Kp-rho*f/sqrt(Kp)*sqrt(vx^2+vr^2))
    src_lin_vx = -visc_int/permeability-rho_int*inertia_coefficient / sqrt(permeability) .* sqrt(vel.vx.^2+vr_int.^2);
end
function src_lin_vr = linear_source_vr(vel, permeability, inertia_coefficient)    
    % Interpolating the radial velocity components, the density and the 
    % viscosity to the position of the axial velocity components
    vx_int = Utilities.interpolate_center2faces((vel.vx(1:end-1,:)+vel.vx(2:end,:))/2,vel.grid.dr,2);
    rho_int = Utilities.interpolate_center2faces(vel.rho,vel.grid.dr,2);
    visc_int = Utilities.interpolate_center2faces(vel.visc,vel.grid.dr,2);

    % Source term = vx*(-μ/Kp-rho*f/sqrt(Kp)*sqrt(vx^2+vr^2))
    src_lin_vr = -visc_int/permeability-rho_int*inertia_coefficient/sqrt(permeability) .* sqrt(vx_int.^2+vel.vr.^2);
end

% Źródła komponentów substancji chemicznych
function src_hydrogen = hydrogen_source(h2, reactions)
    % g/(s*m^3) to kg/(s*m^3)
    unit_conversion = 1/1000;
    src_hydrogen = unit_conversion * h2.M * (3*reactions.R_st+reactions.R_sh);
end
function src_monoxide = monoxide_source(co, reactions)
    % g/(s*m^3) to kg/(s*m^3)
    unit_conversion = 1/1000;
    src_monoxide = unit_conversion * co.M * (reactions.R_st-reactions.R_sh);
end
function src_dioxide = dioxide_source(co2, reactions)
    % g/(s*m^3) to kg/(s*m^3)
    unit_conversion = 1/1000;
    src_dioxide = unit_conversion * co2.M * (reactions.R_sh);
end
function src_methane = methane_source(ch4, reactions)
    % g/(s*m^3) to kg/(s*m^3)
    unit_conversion = 1/1000;
    src_methane = unit_conversion * ch4.M * (-reactions.R_st);
end
function src_water = water_source(h2o, reactions)
    % g/(s*m^3) to kg/(s*m^3)
    unit_conversion = 1/1000;
    src_water = unit_conversion * h2o.M * (-reactions.R_st-reactions.R_sh);
end

% Źródło komponentu energii
function src_q = q_source(temp, reactions, delta_H_st, delta_H_sh)
    src_q = -delta_H_st*reactions.R_st - delta_H_sh*reactions.R_sh;
end

%% Definicja funkcji wyznaczających parametry termofizyczne
% Gęstość - w ChemicalReactions
% Lepkość - w ChemicalReactions
% Pojemność cieplna - w ChemicalReactions
% Dyfuzyjność
function D_eff = effective_species_diffusivity(reactions, species_name, vel_component, temp_component)
    D_eff = (1-sqrt(1-vel_component.porosity)) .* reactions.calculate_species_diffusivity(species_name,vel_component,temp_component);
end
% Przewodność cieplna
function k_eff = effective_thermal_cond(temp_component, reactions, vel_component, solid_phase_thermal_conductivity)
    k_mix = reactions.calculate_mixture_conductivity(temp_component);
    
    k_eff = vel_component.porosity .* k_mix + (1 - vel_component.porosity) .* solid_phase_thermal_conductivity;
end

%% Wprowadzanie członów źródłowych i funkcji wyznaczających parametry termofizyczne

% TODO: zmienić funkcje by korzystały z dependencies, 
% zmienić funkcje by aktualizowały parametry po wprowadzeniu wszystkich zalezności

% Człony źródłowe
vel_component.src_linx_f = @(x) linear_source_vx(x, permeability, inertia_coefficient);
vel_component.src_linr_f = @(x) linear_source_vr(x, permeability, inertia_coefficient);

temp_component.q_function = @(x) q_source(x, reactions, delta_H_st, delta_H_sh);

hydrogen.src_f = @(x) hydrogen_source(x, reactions);
carbon_monoxide.src_f = @(x) monoxide_source(x, reactions);
carbon_dioxide.src_f = @(x) dioxide_source(x, reactions);
methane.src_f = @(x) methane_source(x, reactions);
water.src_f = @(x) water_source(x, reactions);


% Funkcje wyznaczające parametry
vel_component.rho_f = @(x) reactions.calculate_mixture_density(x, temp_component);
vel_component.visc_f = @(x) reactions.calculate_mixture_viscosity(temp_component);

temp_component.cp_function = @(x) reactions.calculate_mixture_heat_capacity(x);
temp_component.k_function = @(x) effective_thermal_cond(x, reactions, vel_component, solid_phase_thermal_conductivity);

hydrogen.D_f = @(x) effective_species_diffusivity(reactions, "hydrogen", vel_component, temp_component);
carbon_monoxide.D_f = @(x) effective_species_diffusivity(reactions, "carbon monoxide", vel_component, temp_component);
carbon_dioxide.D_f = @(x) effective_species_diffusivity(reactions, "carbon dioxide", vel_component, temp_component);
methane.D_f = @(x) effective_species_diffusivity(reactions, "methane", vel_component, temp_component);
water.D_f = @(x) effective_species_diffusivity(reactions, "water", vel_component, temp_component);

%% Definicja warunków brzegowych
% Wartości skalarne

% Wlot
for i=1:sz(2)
    % Określone wartości na wlocie
    temp_component.set_boundary_conditions("west",i,[1 temp_inlet 0]);
    vel_component.set_boundary_conditions("p","west",i,[1 p_atm 0]);
    
    hydrogen.set_boundary_conditions("west",i,[1 0 0]);
    carbon_monoxide.set_boundary_conditions("west",i,[1 0 0]);
    carbon_dioxide.set_boundary_conditions("west",i,[1 0 0]);

    methane.set_boundary_conditions("west",i,[1 M_CH4/(M_CH4+SC*M_H2O) 0]);
    water.set_boundary_conditions("west",i,[1 SC*M_H2O/(M_CH4+SC*M_H2O) 0]);
end
% Wylot
for i=1:sz(2)
    % Brak gradientu na wylocie
    temp_component.set_boundary_conditions("east",i,[2 0 0]);
    vel_component.set_boundary_conditions("p","east",i,[2 0 0]);
    
    hydrogen.set_boundary_conditions("east",i,[2 0 0]);
    carbon_monoxide.set_boundary_conditions("east",i,[2 0 0]);
    carbon_dioxide.set_boundary_conditions("east",i,[2 0 0]);
    methane.set_boundary_conditions("east",i,[2 0 0]);
    water.set_boundary_conditions("east",i,[2 0 0]);
end
% Ściana
for i=1:sz(1)
    temp_component.set_boundary_conditions("north",i,[1 temp_wall 0]); % Stała temperatura ściany
    vel_component.set_boundary_conditions("p","north",i,[2 0 0]); % Brak gradientu na ścianie
    
    % Brak gradientu na ścianie
    hydrogen.set_boundary_conditions("north",i,[2 0 0]);
    carbon_monoxide.set_boundary_conditions("north",i,[2 0 0]);
    carbon_dioxide.set_boundary_conditions("north",i,[2 0 0]);
    methane.set_boundary_conditions("north",i,[2 0 0]);
    water.set_boundary_conditions("north",i,[2 0 0]);
end
% Oś symetrii
for i=1:sz(1)
    % Brak gradientu przy osii symetrii
    temp_component.set_boundary_conditions("south",i,[2 0 0]);
    vel_component.set_boundary_conditions("p","south",i,[2 0 0]);
    
    hydrogen.set_boundary_conditions("south",i,[2 0 0]);
    carbon_monoxide.set_boundary_conditions("south",i,[2 0 0]);
    carbon_dioxide.set_boundary_conditions("south",i,[2 0 0]);
    methane.set_boundary_conditions("south",i,[2 0 0]);
    water.set_boundary_conditions("south",i,[2 0 0]);
end



% Wartości przesunięte w osii X

% Wlot
for i=1:sz(2)
    % Określone wartości na wlocie
    vel_component.set_boundary_conditions("vx","west",i,[1 v_inlet 0]);
end
% Wylot
for i=1:sz(2)
    % Brak gradientu na wylocie
    vel_component.set_boundary_conditions("vx","east",i,[2 0 0]);
end
% Ściana
for i=1:(sz(1)+1)
    % Brak poślizgu na ścianie
    vel_component.set_boundary_conditions("vx","north",i,[1 0 0]);
end
% Oś symetrii
for i=1:(sz(1)+1)
    % Brak gradientu przy osii symetrii
    vel_component.set_boundary_conditions("vx","south",i,[2 0 0]);
end



% Wartości przesunięte w osii R

% Wlot
for i=1:(sz(2)+1)
    % Określone wartości na wlocie
    vel_component.set_boundary_conditions("vr","west",i,[1 0 0]);
end
% Wylot
for i=1:(sz(2)+1)
    % Brak gradientu na wylocie
    vel_component.set_boundary_conditions("vr","east",i,[2 0 0]);
end
% Ściana
for i=1:sz(1)
    % Brak poślizgu na ścianie
    vel_component.set_boundary_conditions("vr","north",i,[1 0 0]);
end
% Oś symetrii
for i=1:sz(1)
    % Brak gradientu przy osii symetrii
    vel_component.set_boundary_conditions("vr","south",i,[2 0 0]);
end

%% Wyznaczanie rozwiązania

% Parametry trzeba najpierw obliczyć kilkukrotnie by otrzymać poprawne wartości
for i = 1:10
    vel_component.update_properties();
    temp_component.update_properties();

    hydrogen.update_properties();
    carbon_monoxide.update_properties();
    carbon_dioxide.update_properties();
    methane.update_properties();
    water.update_properties();

    reactions.update_properties();
end


a = false;

reactions.iterate();

debug_val = false;
while ~a
    sim.noi = sim.noi + 1;

    a = true;

    % Iterowanie każdego komponentu
    a = vel_component.iterate() && a;
    a = temp_component.iterate() && a;

    if debug_val

        a = hydrogen.iterate() && a;
        a = carbon_monoxide.iterate() && a;
        a = carbon_dioxide.iterate() && a;
        a = methane.iterate() && a; 
        a = water.iterate() && a;
    
        a = reactions.iterate() && a;

    end

    % Wyświetlanie rezyduów
    if mod(sim.noi,10)==0 || sim.noi == 1 || true
        fprintf('Current iteration: %d \n',sim.noi);
        fprintf('The vx residual is equal to: %d \n', vel_component.residual_history_vx(sim.noi));
        fprintf('The vr residual is equal to: %d \n', vel_component.residual_history_vy(sim.noi));
        fprintf('The p residual is equal to: %d \n', vel_component.residual_history_p(sim.noi));
        fprintf('The T residual is equal to: %d \n', temp_component.residual_history(sim.noi));

        % fprintf('The H2 residual is equal to: %d \n', hydrogen.residual_history(sim.noi));
        % fprintf('The CO residual is equal to: %d \n', carbon_monoxide.residual_history(sim.noi));
        % fprintf('The CO2 residual is equal to: %d \n', carbon_dioxide.residual_history(sim.noi));
        % fprintf('The CH4 residual is equal to: %d \n', methane.residual_history(sim.noi));
        % fprintf('The H2O residual is equal to: %d \n\n', water.residual_history(sim.noi));
    end
    
    % Zatrzymywanie pętli, jeśli liczba iteracji przekracza dozwoloną wartość
    if sim.noi == max_iters
        break
    end
end

%% Pokazywanie rozwiąń

% Grid
[X,Y]=meshgrid(grid.pos_cent_x,grid.pos_cent_r);

% Interpolated velocities
vx_int=(vel_component.vx(1:end-1,:)+vel_component.vx(2:end,:))/2;
vr_int=(vel_component.vr(:,1:end-1)+vel_component.vr(:,2:end))/2;
v_magn_int=sqrt(vx_int.^2+vr_int.^2);


% Values have to be transposed - probably because MATLAB uses column major order
% Temperature graph
figure(1);
colormap(jet);
contourf(X,Y,transpose(temp_component.temp),30);

% Velocity magnitude graph
figure(2);
colormap(jet);
contourf(X,Y,transpose(v_magn_int),30);

% Velocity in the X direction graph
figure(3);
colormap(jet);
contourf(X,Y,transpose(vx_int),30);

% Velocity in the Y direction graph
figure(4);
colormap(jet);
contourf(X,Y,transpose(vr_int),30);

% Streamlines
figure(5);
streamslice(X,Y,transpose(vx_int),transpose(vr_int),2);

% Residuals
figure(6);
semilogy(vel_component.residual_history_vx(1:sim.noi),'DisplayName','vx residual');
hold on;
semilogy(vel_component.residual_history_vy(1:sim.noi),'DisplayName','vr residual');
semilogy(vel_component.residual_history_p(1:sim.noi),'DisplayName','p residual');
semilogy(temp_component.residual_history(1:sim.noi),'DisplayName','T residual');