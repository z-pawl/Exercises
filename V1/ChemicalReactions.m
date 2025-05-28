classdef ChemicalReactions < IComponent
    properties
        % Dependencies
        sim Simulation2D {mustBeScalarOrEmpty}
        grid Grid2D {mustBeScalarOrEmpty}
        velocity Velocity2DSIMPLE {mustBeScalarOrEmpty}
        temp TemperatureCondConv2DSteady {mustBeScalarOrEmpty}

        % Chemical species
        species dictionary
        species_ind dictionary


        sz double               % Size

        Y_total double          % Sum of mass fractions of all species
        R_st double             % Rate of methane/steam reforming reaction
        R_sh double             % Rate of water-gas-shift reaction

        R double = 8.314        % Universal gas constant

        a double                % Order of reaction with respect to methane
        b double                % Order of reaction with respect to water
        A_st double             % Arrhenius constant
        E_a double              % Activation energy
        delta_G double          % Change of standard Gibbs free energy of water-gas-shift reaction

        SC double               % Steam carbon ratio
        w_cat double            % Catalyst density
        n_inlet_ch4 double      % Volumetric flow rate of methane
        V double                % Reactor volume

        M_mix double            % Molar mass of the mixture

        D_AB double             % Binary diffusion coefficients (without the part depending on temperature and pressure)
        phi_visc double         % Coefficients used to calculate the mixture viscosity
    end
    methods
        function obj = ChemicalReactions(sim, grid, vel, temp, sz,a,b,A_st,E_a,delta_G,SC,w_cat,n_inlet_ch4,V,names, species_components)
            obj.sim = sim;
            obj.grid = grid;
            obj.velocity = vel;
            obj.temp = temp;
            
            obj.sz = sz;

            obj.species = dictionary(names, species_components);
            obj.species_ind = dictionary();

            for i = 1:numel(names)
                obj.species_ind(names(i)) = i;
                species_components(i).reactions = obj;
            end

            obj.a = a;
            obj.b = b;
            obj.A_st = A_st;
            obj.E_a = E_a;
            obj.delta_G = delta_G;

            obj.SC = SC;
            obj.w_cat = w_cat;
            obj.n_inlet_ch4 = n_inlet_ch4;
            obj.V = V;

            obj.M_mix = 1 ./ Utilities.sum_cell(cellfun(@(x) x.Y ./ x.M ,num2cell(values(obj.species)), 'UniformOutput', false));

            obj.calculate_phi_visc();
            obj.calculate_D_AB();
        end
        
        function update_properties(obj)
            
            obj.Y_total = Utilities.sum_cell(cellfun(@(x) x.Y, num2cell(values(obj.species)), 'UniformOutput', false));

            speciess = values(obj.species);
            % Normalizing the mass fractions
            for i = 1:numel(speciess)
                speciess(i).Y = speciess(i).Y ./ obj.Y_total;
            end

            % Chemical reactions rate
            obj.R_st = obj.w_cat * obj.A_st .* exp(-obj.E_a ./ (obj.R .* obj.temp.temp)) ...
                .* obj.species("methane").p_partial.^obj.a .* obj.species("water").p_partial.^obj.b;

            K_sh = exp(-obj.delta_G ./ (obj.R * obj.temp.temp));
            % smol = 1e-50;
            % K_sh = ((obj.species("carbon dioxide").p_partial+smol) .* (obj.species("hydrogen").p_partial+smol)) ./ ((obj.species("carbon monoxide").p_partial+smol) .* (obj.species("water").p_partial+smol));

            xcr = obj.R_st * obj.V / obj.n_inlet_ch4;

            % % Pomysł 1: użyć objętości danej komórki - symulacja nie
            % % dywerguje, jednak również się nie zbiega, Y_CO cały czas jest
            % % równe 0, wskazuje to na to, że w wyniku rozwiązywania równań
            % % Y_CO wychodzi ujemne - zmieniamy to manualnie, powoduje to
            % % fluktuacje i brak zbiegania się pól
            % V_cell = 2*pi*obj.grid.dx.*obj.grid.dr.*obj.grid.pos_cent_r;
            % xcr = obj.R_st .* V_cell / obj.n_inlet_ch4;

            % Coefficients used to calculate ycr
            a_ycr = (K_sh-1);
            b_ycr = -(K_sh*obj.SC+3*xcr);
            c_ycr = K_sh .* xcr .* (obj.SC-xcr);

            ycr = (-b_ycr+sqrt(b_ycr.^2 - 4*a_ycr.*c_ycr)) ./ (2*a_ycr);

            % obj.R_sh = obj.R_st .* ycr;

            % Pomysł 2: R_sh ustawiane na 0 (potem można próbować
            % wyznaczać), a równowaga reakcji zapewniana jest manualnie
            obj.R_sh = zeros(obj.sz);

            % Total molar mass
            obj.M_mix = 1 ./ Utilities.sum_cell(cellfun(@(x) x.Y ./ x.M ,num2cell(values(obj.species)), 'UniformOutput', false));
        end
        
        function converged = iterate(obj)
            obj.update_properties();

            converged=obj.has_converged(0);
        end

        function res = calculate_residual(obj, coeff)
            % Does nothing
            res = 0;
        end

        function converged = has_converged(obj, res)
            converged = true;
        end

        function delete_residual_history(obj)
            % Does nothing
        end

        function phi = phi_k_f(obj, Mi, Mj, visc_i, visc_j)
            M_ratio = Mi/Mj;
            visc_ratio = visc_i / visc_j;

            phi = (1 + (visc_ratio) ^ (0.5) * M_ratio ^ (-0.25)) ^ 2 / (8 * (1 + M_ratio)) ^ (1/2);
        end

        function calculate_D_AB(obj)
            % Constant in the expression for the binary diffusion
            % coefficient
            C = 0.00143;

            species_components = values(obj.species); 
            num_species = numel(species_components);

            for i = 1:num_species
                for j = 1:num_species
                    M_AB = 2/(1/species_components(i).M+1/species_components(j).M);
                    obj.D_AB(i,j) = C / (M_AB^(0.5) ...
                        * ( species_components(i).V_D ^ (1/3) + species_components(j).V_D ^ (1/3))^2);
                end
                % The diagonal is set to infinity
                obj.D_AB(i,i) = double("inf");
            end
        end
        function calculate_phi_visc(obj)
            species_components = values(obj.species); 
            num_species = numel(species_components);

            for i = 1:num_species
                for j = 1:num_species
                    obj.phi_visc(i,j) = (species_components(j).M / species_components(i).M) ^ (0.5);
                end
            end
        end
        
        function D = calculate_species_diffusivity(obj, species_name, vel_component, temp_component)
            ind = obj.species_ind(species_name);

            % Calculating the molar fractions
            molar_fractions = cellfun(@(x) x.Y / x.M .* obj.M_mix, num2cell(values(obj.species)), 'UniformOutput', false);
            
            D = 1 ./ Utilities.sum_cell(cellfun(@(x,y) x ./ y, molar_fractions, num2cell(transpose(obj.D_AB(ind,:))), 'UniformOutput', false));

            % Pascal to bar
            unit_conversion = 1/10^5;
            D = D .* (temp_component.temp .^ 1.75) ./ (vel_component.p * unit_conversion);

            % cm^2/s to m^2/s
            unit_conversion = 1/10000;
            D = D * unit_conversion;
        end
        function rho = calculate_mixture_density(obj, vel_component, temp_component)
            % The density is calculated using the ideal gas law
            rho = vel_component.p .* obj.M_mix ./ (temp_component.temp * obj.R);

            % g/m^3 to kg/m^3
            unit_conversion = 1/1000;
            rho = rho * unit_conversion;
        end
        function cp = calculate_mixture_heat_capacity(obj, temp_component)
            % Calculating the molar fractions
            molar_fractions = cellfun(@(x) x.Y / x.M .* obj.M_mix, num2cell(values(obj.species)),'UniformOutput', false);

            % Calculating the molar heat capacity for each of the species
            molar_cp = cellfun(@(x) arrayfun(@(t) polyval(x.Cp,t), temp_component.temp), num2cell(values(obj.species)), 'UniformOutput', false);

            % Molar heat capacity of the mixture
            cp = Utilities.sum_cell(cellfun(@(x,y) x .* y, molar_cp, molar_fractions, 'UniformOutput',false));

            % Specific heat capacity
            cp = cp .* obj.M_mix;

            % J/(g*K) to J/(kg*K)
            unit_conversion = 1000;
            cp = cp * unit_conversion;
        end
        function visc = calculate_mixture_viscosity(obj, temp_component)
            num_species = numel(values(obj.species));

            % Calculating the molar fractions
            molar_fractions = cellfun(@(x) x.Y / x.M .* obj.M_mix, num2cell(values(obj.species)), 'UniformOutput', false);

            % Calculating the dynamic viscosity for each of the species
            visc_species = cellfun(@(x) arrayfun(@(t) polyval(x.visc,t), temp_component.temp), num2cell(values(obj.species)), 'UniformOutput', false);


            visc = zeros(obj.sz);

            for i = 1:obj.sz(1)
                for j = 1:obj.sz(2)
                    for p = 1:num_species
                        denominator = 0;
                        for q = 1:num_species
                            denominator = denominator + molar_fractions{q}(i,j) * obj.phi_visc(p,q);
                        end
                        visc(i,j) = visc(i,j) + (molar_fractions{p}(i,j) * visc_species{p}(i,j)) / denominator;
                    end
                end
            end
        end
        function k = calculate_mixture_conductivity(obj, temp_component)
            % Used species
            species_comp = num2cell(values(obj.species));

            num_species = numel(species_comp);

            % Calculating the molar fractions
            molar_fractions = cellfun(@(x) x.Y / x.M .* obj.M_mix, species_comp, 'UniformOutput', false);

            % Calculating the dynamic viscosity for each of the species
            visc_species = cellfun(@(x) arrayfun(@(t) polyval(x.visc,t), temp_component.temp), species_comp, 'UniformOutput', false);
            
            % Calculating the thermal conductivity for each of the species
            k_species = cellfun(@(x) arrayfun(@(t) polyval(x.k,t), temp_component.temp), species_comp, 'UniformOutput', false);
        
            k = zeros(obj.sz);

            for i = 1:obj.sz(1)
                for j = 1:obj.sz(2)
                    for p = 1:num_species
                        denominator = 0;
                        for q = 1:num_species
                            denominator = denominator + ...
                            molar_fractions{q}(i,j) * ...
                            obj.phi_k_f(species_comp{p}.M, species_comp{q}.M, visc_species{p}(i,j), visc_species{q}(i,j));
                        end
                        k(i,j) = k(i,j) + molar_fractions{p}(i,j) * k_species{p}(i,j) / denominator;
                    end
                end
            end
        end
    end
end