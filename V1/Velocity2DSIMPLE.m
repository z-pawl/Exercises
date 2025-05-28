% Pomysły na miejsce umieszczenia podrelaksacji
% -przy przypisywaniu nowej prędkości - trochę wolniejsze, mniejszy wsp
% relaksacji ciśnienia (rzędu 0.1)
% -przy rozwiązywaniu równań - szybsze, większy wsp relaksacji ciśnienia
% Porównać stabilność

% Opłaca się bardziej rozwiązywać równanie ciśnienia - może multigrid
% będzie dobrym rozwiązaniem

classdef Velocity2DSIMPLE < IComponent
    properties
        % Simulation object
        sim Simulation2D {mustBeScalarOrEmpty}

        % Dependencies
        grid Grid2D {mustBeScalarOrEmpty}

        sz double % Size

        % A staggered grid approach is used to represent the velocity field
        vx double
        vr double
        % A collocated grid approach is used to represent the pressure field
        p double
        
        porosity double      % Porosity
        rho double           % Density
        visc double          % Dynamic viscosity
        srcx double          % Constant source term in the x direction
        src_linx double      % Linear source term in the x direction S_total = S_const + S_linear * vx
        srcr double         % Constant source term in the r direction
        src_linr double     % Linear source term in the r direction S_total = S_const + S_linear * vr

        % Functions used to calculate the physical properties
        rho_f (1,1) function_handle = @(x) 0
        visc_f (1,1) function_handle = @(x) 0
        srcx_f (1,1) function_handle = @(x) 0
        src_linx_f (1,1) function_handle = @(x) 0
        srcr_f (1,1) function_handle = @(x) 0
        src_linr_f (1,1) function_handle = @(x) 0

        tol_vel (1,1) double {mustBePositive} = 1e-7
        tol_pressure (1,1)  double {mustBePositive} = 1e-7
        residual_history_vx (:,1) double = zeros(10000,1)
        residual_history_vy (:,1) double = zeros(10000,1)
        residual_history_p (:,1) double = zeros(10000,1)
        rel_fact_v (1,1) double = 1
        rel_fact_p (1,1) double = 0.8

        % boundaries_vx IBoundary
        % boundaries_vy IBoundary
        % boundaries_p IBoundary

        boundaries_west_vx (:,3)
        boundaries_east_vx (:,3)
        boundaries_south_vx (:,3)
        boundaries_north_vx (:,3)

        boundaries_west_vr (:,3)
        boundaries_east_vr (:,3)
        boundaries_south_vr (:,3)
        boundaries_north_vr (:,3)

        boundaries_west_p (:,3)
        boundaries_east_p (:,3)
        boundaries_south_p (:,3)
        boundaries_north_p (:,3)

    end
    methods
        function obj = Velocity2DSIMPLE(sim, grid, sz, vx, vr, p, porosity, rho_f, visc_f, srcx_f, src_linx_f, srcr_f, src_linr_f, tol_vel, tol_pressure, rel_fact_v, rel_fact_p)
            
            if ~isequal(sz+[1 0], size(vx)) || ~isequal(sz+[0 1], size(vr)) || ~isequal(sz, size(p))
                error("Arrays provided in the input must have the size conforming with the specified size");
            end

            obj.sim = sim;
            obj.grid = grid;
            
            obj.sz=sz;

            obj.vx=vx;
            obj.vr=vr;
            obj.p=p;

            obj.porosity = porosity;

            obj.rho_f=rho_f;
            obj.visc_f=visc_f;
            obj.srcx_f=srcx_f;
            obj.src_linx_f=src_linx_f;
            obj.srcr_f=srcr_f;
            obj.src_linr_f=src_linr_f;

            obj.tol_vel = tol_vel;
            obj.tol_pressure = tol_pressure;
            obj.rel_fact_v = rel_fact_v;
            obj.rel_fact_p = rel_fact_p;

            obj.boundaries_west_vx=zeros(obj.sz(2),3);
            obj.boundaries_east_vx=zeros(obj.sz(2),3);
            obj.boundaries_south_vx=zeros(obj.sz(1)+1,3);
            obj.boundaries_north_vx=zeros(obj.sz(1)+1,3);

            obj.boundaries_west_vr=zeros(obj.sz(2)+1,3);
            obj.boundaries_east_vr=zeros(obj.sz(2)+1,3);
            obj.boundaries_south_vr=zeros(obj.sz(1),3);
            obj.boundaries_north_vr=zeros(obj.sz(1),3);

            obj.boundaries_west_p=zeros(obj.sz(2),3);
            obj.boundaries_east_p=zeros(obj.sz(2),3);
            obj.boundaries_south_p=zeros(obj.sz(1),3);
            obj.boundaries_north_p=zeros(obj.sz(1),3);
        end

        function update_properties(obj)
            obj.rho=obj.rho_f(obj);
            obj.visc=obj.visc_f(obj);
            obj.srcx=obj.srcx_f(obj);
            obj.src_linx=obj.src_linx_f(obj);
            obj.srcr=obj.srcr_f(obj);
            obj.src_linr=obj.src_linr_f(obj);
        end

        function [coeff_vx, coeff_vr] = get_coefficients_vel(obj)
            coeff_vx=obj.get_coefficients_vx();
            coeff_vr=obj.get_coefficients_vr();

            % Pressure force is already included in the source term
            [A_dpx, A_dpr] = obj.get_pressure_force(obj.p,false);

            coeff_vx(:,:,6)=coeff_vx(:,:,6)+A_dpx;
            coeff_vr(:,:,6)=coeff_vr(:,:,6)+A_dpr;
        end

        function coeff = get_coefficients_vx(obj)
            % Preallocating the coefficient matrix
            coeff = zeros([(obj.sz+[1 0]) 6]);

            % A row of zeros used to pad the data in order to make sure that there is no index out of bounds error
            row0 = zeros(1,obj.sz(2));

            area_x = 2*pi*obj.grid.dr.*obj.grid.pos_cent_r;
            area_r = 2*pi*obj.grid.dx_stag.*obj.grid.pos_stag_r;
            volume = 2*pi*obj.grid.dx_stag.*obj.grid.dr.*obj.grid.pos_cent_r;

            % Vectors holding values of F, D and A(|P|) for every interface
            % between control volumes along x and r axis respectively
            % (for example Fx(1,1) is a value for the interface between the control volumes with indexes (1,1) and (2,1) respectively)

            % Flow rate through the cell face: F = rho*v*A
            % The density is calculated as the mean of the densities of the neighbouring cells
            % The velocity is taken as the velocity at the cell face
            % A - area of the cell face

            % The expression for Fx at the eastern and western boundary does not exist as there are no scalar value cells on the other side of the boundary
            % There is no need to interpolate the density as the staggered cell cell face center coincides with the scalar value cell center
            % Nor is there any need to lineary interpolate the velocity as the distance from the staggered cell face to the staggered cell center on either side of the face is the same
            Fx=[obj.rho(1,:).*obj.vx(1,:); obj.rho.*(obj.vx(1:end-1,:)+obj.vx(2:end,:))/2; obj.rho(end,:).*obj.vx(end,:)] .* area_x;

            % Densities interpolated to the cell faces along the r direction
            rho_int_r = Utilities.interpolate_center2faces(Utilities.interpolate_center2faces(obj.rho,obj.grid.dr,2),obj.grid.dx,1);
            
            % Velocities in the r direction interpolated to the cell faces along the r direction
            vr_int_r = Utilities.interpolate_center2faces(obj.vr,obj.grid.dx,1);

            Fr = rho_int_r .* vr_int_r .* area_r;

            % Conductance: D =  k/dx*A

            % There is no need to interpolate the viscosity as the staggered cell cell face center coincides with the scalar value cell center
            Dx=[row0; obj.visc./obj.grid.dx; row0] .* area_x;

            Dr = Utilities.interpolate_center2faces(Utilities.interpolate_center2faces(obj.visc,obj.grid.dr,2),obj.grid.dx,1) ...
                ./ obj.grid.dr_stag .* area_r;

            
            % Taking porosity into account
            Fx = Fx / obj.porosity^2;
            Fr = Fr / obj.porosity^2;
            Dx = Dx / obj.porosity;
            Dr = Dr / obj.porosity;
            

            % Function A(|P|) calculated using the power law scheme
            Ax = [row0; (1-0.1*abs(Fx(2:end-1,:)./Dx(2:end-1,:))).^5; row0];
            Ar = (1-0.1*abs(Fr./Dr)).^5;

            % Every negative value is changed to 0
            Ax = Ax.*(Ax>0);
            Ar = Ar.*(Ar>0);

            % aW = Dw*A(|Pw|)+max(Fw,0)
            coeff(:,:,1) = Dx(1:end-1,:).*Ax(1:end-1,:)+Fx(1:end-1,:).*(Fx(1:end-1,:)>0);
            % aE = De*A(|Pe|)+max(-Fe,0)
            coeff(:,:,2) = Dx(2:end,:).*Ax(2:end,:)-Fx(2:end,:).*(Fx(2:end,:)<0);
            % aS = Ds*A(|Ps|)+max(Fs,0)
            coeff(:,:,3) = Dr(:,1:end-1).*Ar(:,1:end-1)+Fr(:,1:end-1).*(Fr(:,1:end-1)>0);
            % aN = Dn*A(|Pn|)+max(-Fn,0)
            coeff(:,:,4) = Dr(:,2:end).*Ar(:,2:end)-Fr(:,2:end).*(Fr(:,2:end)<0);
            % aP = aE+aW+aN+aS-SP*dx*dy
            coeff(:,:,5) = coeff(:,:,1)+coeff(:,:,2)+coeff(:,:,3)+coeff(:,:,4)-obj.src_linx(:,:).*volume+(Fx(2:end,:)-Fx(1:end-1,:))+(Fr(:,2:end)-Fr(:,1:end-1));
            % b = SC*dx*dy
            coeff(:,:,6) = obj.srcx(:,:).*volume;

            coeff = obj.apply_boundary_conditions_vx(coeff);
        end
        function coeff = get_coefficients_vr(obj)
            % Preallocating the coefficient matrix
            coeff = zeros([(obj.sz+[0 1]) 6]);

            % A column of zeros used to pad the data in order to make sure that there is no index out of bounds error
            col0 = zeros(obj.sz(1),1);

            area_x = 2*pi*obj.grid.dr_stag.*([obj.grid.pos_stag_r(1,1) obj.grid.pos_cent_r] + [obj.grid.pos_cent_r obj.grid.pos_stag_r(1,end)])/2;
            area_r = 2*pi*obj.grid.dx.*[obj.grid.pos_stag_r(1,1) obj.grid.pos_cent_r obj.grid.pos_stag_r(1,end)];
            volume = 2*pi*obj.grid.dx.*obj.grid.dr_stag.*([obj.grid.pos_stag_r(1,1) obj.grid.pos_cent_r] + [obj.grid.pos_cent_r obj.grid.pos_stag_r(1,end)])/2;

            % Vectors holding values of F, D and A(|P|) for every interface
            % between control volumes along x and r axis respectively
            % (for example Fx(1,1) is a value for the interface between the control volumes with indexes (1,1) and (2,1) respectively)

            % Flow rate through the cell face: F = rho*v*A
            % The density is calculated as the mean of the densities of the neighbouring cells
            % The velocity is taken as the velocity at the cell face
            % A - area of the cell face

            Fx = Utilities.interpolate_center2faces(Utilities.interpolate_center2faces(obj.rho, obj.grid.dx, 1), obj.grid.dr, 2) ...
                .* Utilities.interpolate_center2faces(obj.vx, obj.grid.dr, 2) .* area_x;

            % There is no need to interpolate the density as the staggered cell cell face center coincides with the scalar value cell center
            % Nor is there any need to lineary interpolate the velocity as the distance from the staggered cell face to the staggered cell center on either side of the face is the same
            Fr=[obj.rho(:,1).*obj.vr(:,1) obj.rho.*(obj.vr(:,1:end-1)+obj.vr(:,2:end))/2 obj.rho(:,end).*obj.vr(:,end)] .* area_r;

            % Conductance: D =  k/dx*A

            Dx=Utilities.interpolate_center2faces(Utilities.interpolate_center2faces(obj.visc, obj.grid.dx, 1), obj.grid.dr, 2) ...
                ./ obj.grid.dx_stag .* area_x;

            % There is no need to interpolate the viscosity as the staggered cell cell face center coincides with the scalar value cell center
            Dr = [col0 obj.visc./obj.grid.dr col0] .* area_r;
            
            % Taking porosity into account
            Fx = Fx / obj.porosity^2;
            Fr = Fr / obj.porosity^2;
            Dx = Dx / obj.porosity;
            Dr = Dr / obj.porosity;

            % Function A(|P|) calculated using the power law scheme
            Ax = (1-0.1*abs(Fx./Dx)).^5;
            Ar = [col0 (1-0.1*abs(Fr(:,2:end-1)./Dr(:,2:end-1))).^5 col0];

            % Every negative value is changed to 0
            Ax = Ax.*(Ax>0);
            Ar = Ar.*(Ar>0);

            % aW = Dw*A(|Pw|)+max(Fw,0)
            coeff(:,:,1) = Dx(1:end-1,:).*Ax(1:end-1,:)+Fx(1:end-1,:).*(Fx(1:end-1,:)>0);
            % aE = De*A(|Pe|)+max(-Fe,0)
            coeff(:,:,2) = Dx(2:end,:).*Ax(2:end,:)-Fx(2:end,:).*(Fx(2:end,:)<0);
            % aS = Ds*A(|Ps|)+max(Fs,0)
            coeff(:,:,3) = Dr(:,1:end-1).*Ar(:,1:end-1)+Fr(:,1:end-1).*(Fr(:,1:end-1)>0);
            % aN = Dn*A(|Pn|)+max(-Fn,0)
            coeff(:,:,4) = Dr(:,2:end).*Ar(:,2:end)-Fr(:,2:end).*(Fr(:,2:end)<0);
            % aP = aE+aW+aN+aS-SP*dx*dy
            coeff(:,:,5) = coeff(:,:,1)+coeff(:,:,2)+coeff(:,:,3)+coeff(:,:,4) ...
                +Utilities.interpolate_center2faces(obj.visc, obj.grid.dr, 2).* volume./(obj.porosity .* obj.grid.pos_stag_r .^ 2) -obj.src_linr(:,:).*volume+(Fx(2:end,:) ...
                -Fx(1:end-1,:))+(Fr(:,2:end)-Fr(:,1:end-1));
            % b = SC*dx*dy
            coeff(:,:,6) = obj.srcr(:,:).*volume;

            coeff = obj.apply_boundary_conditions_vr(coeff);
        end
        function coeff = get_coefficients_p_corr(obj, coeff_vx_aP, coeff_vy_aP, vx_star, vr_star)
            % Preallocating the coefficient matrix
            coeff = zeros([(obj.sz) 6]);

            area_x = 2*pi*obj.grid.dr.*obj.grid.pos_cent_r;
            area_r = 2*pi*obj.grid.dx.*obj.grid.pos_stag_r;

            % Density interpolated to the cell faces along the x direction
            rho_int_x = Utilities.interpolate_center2faces(obj.rho, obj.grid.dx, 1);
                
            % Density interpolated to the cell faces along the r direction
            rho_int_r = Utilities.interpolate_center2faces(obj.rho, obj.grid.dr, 2);

            % aE = rho_e*d_e*A_e=rho_e*A_e^2/(aP_vx_e);
            aX=rho_int_x.*(area_x.^2)./coeff_vx_aP;
            aR=rho_int_r.*(area_r.^2)./coeff_vy_aP;

            % Hence
            % aW
            coeff(:,:,1)=aX(1:end-1,:);
            % aE
            coeff(:,:,2)=aX(2:end,:);
            % aS
            coeff(:,:,3)=aR(:,1:end-1);
            % aN
            coeff(:,:,4)=aR(:,2:end);
            % aP
            coeff(:,:,5)=coeff(:,:,1)+coeff(:,:,2)+coeff(:,:,3)+coeff(:,:,4);

            % b = ((rho*u*A)_w-(rho*u*A)_e)+((rho*v*A)_s-(rho*v*A)_n)
            rho_vx_A=rho_int_x.*vx_star.*area_x;
            rho_vr_A=rho_int_r.*vr_star.*area_r;

            coeff(:,:,6)=(rho_vx_A(1:end-1,:)-rho_vx_A(2:end,:))+(rho_vr_A(:,1:end-1)-rho_vr_A(:,2:end));

            coeff=obj.apply_boundary_conditions_p_corr(coeff);
        end
        function [A_delta_p_x, A_delta_p_r] = get_pressure_force(obj, p, is_pressure_correction)
            area_x = 2*pi*obj.grid.dr.*obj.grid.pos_cent_r;
            area_r = 2*pi*obj.grid.dx.*obj.grid.pos_stag_r;

            % Boundary conditions are used to obtain the values of the pressure at the boundary
            % Pressure is padded with additional zeros to make place for the values at the boundary
            p = obj.apply_boundary_conditions_p_force(paddata(p,size(p)+[2 2],"Side","both"),is_pressure_correction);

            A_delta_p_x = (p(1:end-1,2:end-1)-p(2:end,2:end-1)) .* area_x;
            A_delta_p_r = (p(2:end-1,1:end-1)-p(2:end-1,2:end)) .* area_r;
        end
        
        function coeff = apply_boundary_conditions_vx(obj, coeff)
            % Function calculating the changes in the coefficient for a given boundary of a cell 
            function [d_aP, d_b] = apply_boundary_condtion_cell(bd_data, aB, A, n_dir)
                d_aP=zeros(size(bd_data,1),1);
                d_b=zeros(size(bd_data,1),1);
                for i=1:size(bd_data,1)
                    if n_dir(i,1)~=0
                        % Half-cell boundary
                        if bd_data(i,1)==1
                            % Dirichlet (given temperature) boundary condition
                            d_aP(i)=-aB(i)+1e30;
                            d_b(i)=1e30*bd_data(i,2);
                        elseif bd_data(i,1)==2
                            % Neumann (given heat flux) boundary condition
                            d_aP(i)=-aB(i);
                            d_b(i)=bd_data(i,2)*A(i)*n_dir(i,1);
                        else
                            % DOES NOT WORK
                            % Robin (convective boundary) boundary condition
                            % m=bd_data(i,3)*dq+2*k;
                            % d_aP=-aB*2*k/m;
                            % d_b=aB*(m-2*k)/m*bd_data(2);
                        end
                    else
                        % Normal-cell boundary
                        if bd_data(i,1)==1
                            % Dirichlet (given temperature) boundary condition
                            d_aP(i)=0;
                            d_b(i)=bd_data(i,2)*aB(i);
                        elseif bd_data(i,1)==2
                            % Neumann (given heat flux) boundary condition
                            d_aP(i)=-aB(i);
                            d_b(i)=bd_data(i,2)*A(i)*n_dir(i,2);
                        else
                            % DOES NOT WORK
                            % Robin (convective boundary) boundary condition
                            % m=bd_data(i,3)*dq+2*k;
                            % d_aP=-aB*2*k/m;
                            % d_b=aB*(m-2*k)/m*bd_data(2);
                        end
                    end
                end
            end

            area_x = 2*pi*obj.grid.dr.*obj.grid.pos_cent_r;
            area_r = 2*pi*obj.grid.dx_stag.*obj.grid.pos_stag_r;
            
            % Calculations of the boundary conditions are performed for every wall separately
            [d_aP_west, d_b_west]=apply_boundary_condtion_cell(obj.boundaries_west_vx, coeff(1,:,1), area_x(1,:), [ones(obj.sz(2),1) zeros(obj.sz(2),1)]);
            
            [d_aP_east, d_b_east]=apply_boundary_condtion_cell(obj.boundaries_east_vx, coeff(end,:,2), area_x(end,:), [-ones(obj.sz(2),1) zeros(obj.sz(2),1)]);
            
            [d_aP_south, d_b_south]=apply_boundary_condtion_cell(obj.boundaries_south_vx, coeff(:,1,3), area_r(:,1), [zeros(obj.sz(1)+1,1) ones(obj.sz(1)+1,1)]);

            [d_aP_north, d_b_north]=apply_boundary_condtion_cell(obj.boundaries_north_vx, coeff(:,end,4), area_r(:,end), [zeros(obj.sz(1)+1,1) -ones(obj.sz(1)+1,1)]);

            % Obtained values are used to change the coefficients of the discretized equations
            % West
            coeff(1,:,5)=coeff(1,:,5)+reshape(d_aP_west,[1 obj.sz(2)]);
            coeff(1,:,6)=coeff(1,:,6)+reshape(d_b_west,[1 obj.sz(2)]);
            coeff(1,:,1)=0;
            % East
            coeff(end,:,5)=coeff(end,:,5)+reshape(d_aP_east,[1 obj.sz(2)]);
            coeff(end,:,6)=coeff(end,:,6)+reshape(d_b_east,[1 obj.sz(2)]);
            coeff(end,:,2)=0;
            % South
            coeff(:,1,5)=coeff(:,1,5)+d_aP_south;
            coeff(:,1,6)=coeff(:,1,6)+d_b_south;
            coeff(:,1,3)=0;
            % North
            coeff(:,end,5)=coeff(:,end,5)+d_aP_north;
            coeff(:,end,6)=coeff(:,end,6)+d_b_north;
            coeff(:,end,4)=0;
        end
        function coeff = apply_boundary_conditions_vr(obj, coeff)
            % Function calculating the changes in the coefficient for a given boundary of a cell 
            function [d_aP, d_b] = apply_boundary_condtion_cell(bd_data, aB, A, n_dir)
                d_aP=zeros(size(bd_data,1),1);
                d_b=zeros(size(bd_data,1),1);
                for i=1:size(bd_data,1)
                    if n_dir(i,1)~=0
                        % Normal-cell boundary
                        if bd_data(i,1)==1
                            % Dirichlet (given temperature) boundary condition
                            d_aP(i)=0;
                            d_b(i)=bd_data(i,2)*aB(i);
                        elseif bd_data(i,1)==2
                            % Neumann (given heat flux) boundary condition
                            d_aP(i)=-aB(i);
                            d_b(i)=bd_data(i,2)*A(i)*n_dir(i,1);
                        else
                            % DOES NOT WORK
                            % Robin (convective boundary) boundary condition
                            % m=bd_data(i,3)*dq+2*k;
                            % d_aP=-aB*2*k/m;
                            % d_b=aB*(m-2*k)/m*bd_data(2);
                        end
                    else
                        % Half-cell boundary
                        if bd_data(i,1)==1
                            % Dirichlet (given temperature) boundary condition
                            d_aP(i)=-aB(i)+1e30;
                            d_b(i)=1e30*bd_data(i,2);
                        elseif bd_data(i,1)==2
                            % Neumann (given heat flux) boundary condition
                            d_aP(i)=-aB(i);
                            d_b(i)=bd_data(i,2)*A(i)*n_dir(i,2);
                        else
                            % DOES NOT WORK
                            % Robin (convective boundary) boundary condition
                            % m=bd_data(i,3)*dq+2*k;
                            % d_aP=-aB*2*k/m;
                            % d_b=aB*(m-2*k)/m*bd_data(2);
                        end
                    end
                end
            end

            area_x = 2*pi*obj.grid.dr_stag.*([obj.grid.pos_stag_r(1,1) obj.grid.pos_cent_r] + [obj.grid.pos_cent_r obj.grid.pos_stag_r(1,end)])/2;
            area_r = 2*pi*obj.grid.dx.*[obj.grid.pos_stag_r(1,1) obj.grid.pos_cent_r obj.grid.pos_stag_r(1,end)];
            
            % Calculations of the boundary conditions are performed for every wall separately
            [d_aP_west, d_b_west]=apply_boundary_condtion_cell(obj.boundaries_west_vr, coeff(1,:,1), area_x(1,:), [ones(obj.sz(2)+1,1) zeros(obj.sz(2)+1,1)]);
            
            [d_aP_east, d_b_east]=apply_boundary_condtion_cell(obj.boundaries_east_vr, coeff(end,:,2), area_x(end,:), [-ones(obj.sz(2)+1,1) zeros(obj.sz(2)+1,1)]);
            
            [d_aP_south, d_b_south]=apply_boundary_condtion_cell(obj.boundaries_south_vr, coeff(:,1,3), area_r(:,1), [zeros(obj.sz(1),1) ones(obj.sz(1),1)]);

            [d_aP_north, d_b_north]=apply_boundary_condtion_cell(obj.boundaries_north_vr, coeff(:,end,4), area_r(:,end), [zeros(obj.sz(1),1) -ones(obj.sz(1),1)]);

            % Obtained values are used to change the coefficients of the discretized equations
            % West
            coeff(1,:,5)=coeff(1,:,5)+reshape(d_aP_west, [1, obj.sz(2)+1]);
            coeff(1,:,6)=coeff(1,:,6)+reshape(d_b_west, [1, obj.sz(2)+1]);
            coeff(1,:,1)=0;
            % East
            coeff(end,:,5)=coeff(end,:,5)+reshape(d_aP_east, [1, obj.sz(2)+1]);
            coeff(end,:,6)=coeff(end,:,6)+reshape(d_b_east, [1, obj.sz(2)+1]);
            coeff(end,:,2)=0;
            % South
            coeff(:,1,5)=coeff(:,1,5)+d_aP_south;
            coeff(:,1,6)=coeff(:,1,6)+d_b_south;
            coeff(:,1,3)=0;
            % North
            coeff(:,end,5)=coeff(:,end,5)+d_aP_north;
            coeff(:,end,6)=coeff(:,end,6)+d_b_north;
            coeff(:,end,4)=0;
        end
        function coeff = apply_boundary_conditions_p_corr(obj, coeff)
            % Function calculating the changes in the coefficient for a given boundary of a cell 
            function [d_aP, d_b] = apply_boundary_condtion_cell(bd_data, aB)
                d_aP=zeros(size(bd_data,1),1);
                d_b=zeros(size(bd_data,1),1);
                for i=1:size(bd_data,1)
                    % Normal-cell boundary
                    if bd_data(i,1)==1
                        % Dirichlet (given temperature) boundary condition
                        d_aP(i)=0;
                        d_b(i)=0;
                    elseif bd_data(i,1)==2
                        % Neumann (given heat flux) boundary condition
                        d_aP(i)=-aB(i);
                        d_b(i)=0;
                    else
                        % DOES NOT WORK
                        % Robin (convective boundary) boundary condition
                        % m=bd_data(i,3)*dq+2*k;
                        % d_aP=-aB*2*k/m;
                        % d_b=aB*(m-2*k)/m*bd_data(2);
                    end
                end
            end
            
            % Calculations of the boundary conditions are performed for every wall separately
            [d_aP_west, d_b_west]=apply_boundary_condtion_cell(obj.boundaries_west_p, coeff(1,:,1));
            
            [d_aP_east, d_b_east]=apply_boundary_condtion_cell(obj.boundaries_east_p, coeff(end,:,2));
            
            [d_aP_south, d_b_south]=apply_boundary_condtion_cell(obj.boundaries_south_p, coeff(:,1,3));

            [d_aP_north, d_b_north]=apply_boundary_condtion_cell(obj.boundaries_north_p, coeff(:,end,4));

            % Obtained values are used to change the coefficients of the discretized equations
            % West
            coeff(1,:,5)=coeff(1,:,5)+reshape(d_aP_west, [1, obj.sz(2)]);
            coeff(1,:,6)=coeff(1,:,6)+reshape(d_b_west, [1, obj.sz(2)]);
            coeff(1,:,1)=0;
            % East
            coeff(end,:,5)=coeff(end,:,5)+reshape(d_aP_east, [1, obj.sz(2)]);
            coeff(end,:,6)=coeff(end,:,6)+reshape(d_b_east, [1, obj.sz(2)]);
            coeff(end,:,2)=0;
            % South
            coeff(:,1,5)=coeff(:,1,5)+d_aP_south;
            coeff(:,1,6)=coeff(:,1,6)+d_b_south;
            coeff(:,1,3)=0;
            % North
            coeff(:,end,5)=coeff(:,end,5)+d_aP_north;
            coeff(:,end,6)=coeff(:,end,6)+d_b_north;
            coeff(:,end,4)=0;
        end
        function p = apply_boundary_conditions_p_force(obj, p, is_pressure_correction)
            function p_b = apply_boundary_condition_cell(bd_data,p_P,n_dir,dq)
                if is_pressure_correction
                    switch bd_data(1)
                        case 1
                            p_b=0;
                        case 2
                            p_b=p_P;
                        case 3
                            % DOES NOT WORK
                    end
                else
                    switch bd_data(1)
                        case 1
                            p_b=bd_data(2);
                        case 2
                            p_b=p_P-bd_data(2)*dq/(2*n_dir);
                        case 3
                            % DOES NOT WORK
                    end
                end
            end

            % Calculations of the boundary conditions are performed for every wall separately
            p_b_west=arrayfun(@(x) apply_boundary_condition_cell(obj.boundaries_west_p(x,:),p(2,x+1),1,obj.grid.dx(1)), 1:obj.sz(2));
            
            p_b_east=arrayfun(@(x) apply_boundary_condition_cell(obj.boundaries_east_p(x,:),p(end-1,x+1),-1,obj.grid.dx(end)), 1:obj.sz(2));
            
            p_b_south=arrayfun(@(x) apply_boundary_condition_cell(obj.boundaries_south_p(x,:),p(x+1,2),1,obj.grid.dr(1)), 1:obj.sz(1));

            p_b_north=arrayfun(@(x) apply_boundary_condition_cell(obj.boundaries_north_p(x,:),p(x+1,end-1),-1,obj.grid.dr(end)), 1:obj.sz(1));
        
            
            % West
            p(1,2:end-1)=p_b_west;
            % East
            p(end,2:end-1)=p_b_east;
            % South
            p(2:end-1,1)=reshape(p_b_south,[obj.sz(1) 1]);
            % North
            p(2:end-1,end)=reshape(p_b_north, [obj.sz(1),1]);
        end
        
        function obj = set_boundary_conditions(obj, field_name, boundary_side, boundary_index, bd_data)
            % Checking the validity of the boundary data
            if ~ismember(bd_data(1), [1 2 3])
                disp("First value of the boundary data describes the boundary type; it has to be either 1, 2 or 3");
                return;
            end

            switch field_name
                case "vx"
                    % Changing the corresponding matrix holding the boundary condition
                    switch boundary_side
                        case "west"
                            if 1<=boundary_index && boundary_index<=obj.sz(2)
                                obj.boundaries_west_vx(boundary_index,:)=bd_data;
                            else
                                disp("Boundary index is out of range");
                            end
                        case "east"
                            if 1<=boundary_index && boundary_index<=obj.sz(2)
                                obj.boundaries_east_vx(boundary_index,:)=bd_data;
                            else
                                disp("Boundary index is out of range");
                            end
                        case "south"
                            if 1<=boundary_index && boundary_index<=obj.sz(1)+1
                                obj.boundaries_south_vx(boundary_index,:)=bd_data;
                            else
                                disp("Boundary index is out of range");
                            end
                        case "north"
                            if 1<=boundary_index && boundary_index<=obj.sz(1)+1
                                obj.boundaries_north_vx(boundary_index,:)=bd_data;
                            else
                                disp("Boundary index is out of range");
                            end
                        otherwise
                            disp("Invalid boundary side name");
                    end
                case "vr"
                    % Changing the corresponding matrix holding the boundary condition
                    switch boundary_side
                        case "west"
                            if 1<=boundary_index && boundary_index<=obj.sz(2)+1
                                obj.boundaries_west_vr(boundary_index,:)=bd_data;
                            else
                                disp("Boundary index is out of range");
                            end
                        case "east"
                            if 1<=boundary_index && boundary_index<=obj.sz(2)+1
                                obj.boundaries_east_vr(boundary_index,:)=bd_data;
                            else
                                disp("Boundary index is out of range");
                            end
                        case "south"
                            if 1<=boundary_index && boundary_index<=obj.sz(1)
                                obj.boundaries_south_vr(boundary_index,:)=bd_data;
                            else
                                disp("Boundary index is out of range");
                            end
                        case "north"
                            if 1<=boundary_index && boundary_index<=obj.sz(1)
                                obj.boundaries_north_vr(boundary_index,:)=bd_data;
                            else
                                disp("Boundary index is out of range");
                            end
                        otherwise
                            disp("Invalid boundary side name");
                    end
                case "p"
                    % Changing the corresponding matrix holding the boundary condition
                    switch boundary_side
                        case "west"
                            if 1<=boundary_index && boundary_index<=obj.sz(2)
                                obj.boundaries_west_p(boundary_index,:)=bd_data;
                            else
                                disp("Boundary index is out of range");
                            end
                        case "east"
                            if 1<=boundary_index && boundary_index<=obj.sz(2)
                                obj.boundaries_east_p(boundary_index,:)=bd_data;
                            else
                                disp("Boundary index is out of range");
                            end
                        case "south"
                            if 1<=boundary_index && boundary_index<=obj.sz(1)
                                obj.boundaries_south_p(boundary_index,:)=bd_data;
                            else
                                disp("Boundary index is out of range");
                            end
                        case "north"
                            if 1<=boundary_index && boundary_index<=obj.sz(1)
                                obj.boundaries_north_p(boundary_index,:)=bd_data;
                            else
                                disp("Boundary index is out of range");
                            end
                        otherwise
                            disp("Invalid boundary side name");
                    end
            end
        end

        function converged = iterate(obj)
            % Updating the properties
            obj.update_properties();


            % Calculating the coefficients of momentum discretized equations and solving them
            [coeff_vx, coeff_vr]=obj.get_coefficients_vel();

            % A variable used for debugging to check if any of the central cell coefficients is equal to zero
            % which would result in a division by zero when solving the equations
            x = any(any(coeff_vx(:,:,5)==0)) || any(any(coeff_vr(:,:,5)==0));

            % EXPLICIT RELAXATION
            % Starred velocities
            vx_star=solve(coeff_vx,obj.vx,1,[randi([1 4])],obj.rel_fact_v);
            vr_star=solve(coeff_vr,obj.vr,1,[randi([1 4])],obj.rel_fact_v);


            % Calculating the coefficients of the pressure correction equation and solving it
            coeff_p_corr=obj.get_coefficients_p_corr(coeff_vx(:,:,5),coeff_vr(:,:,5),vx_star,vr_star);

            p_corr=solve(coeff_p_corr,zeros(obj.sz),50,[1;2;3;4],1);

            % Calculating the pressure force from the pressure correction and using it to correct the velocities
            [A_delta_p_x, A_delta_p_r]=obj.get_pressure_force(p_corr, true);

            vx_star=vx_star+A_delta_p_x./coeff_vx(:,:,5);
            vr_star=vr_star+A_delta_p_r./coeff_vr(:,:,5);

            % Updating the fields using the relaxation factor
            obj.p=obj.p+p_corr*obj.rel_fact_p;
            obj.vx=vx_star;
            obj.vr=vr_star;

            % Updating properties and calculating the coefficients of the discretized equations used to calculate residuals
            obj.update_properties();
            [coeff_vx, coeff_vr]=obj.get_coefficients_vel();
            coeff_p_corr=obj.get_coefficients_p_corr(coeff_vx(:,:,5),coeff_vr(:,:,5),obj.vx,obj.vr);

            % Calculating residual and saving them
            res_vx=obj.calculate_residual(coeff_vx,"vx");
            res_vr=obj.calculate_residual(coeff_vr,"vr");
            res_p=obj.calculate_residual(coeff_p_corr,"p");
            
            obj.residual_history_vx(obj.sim.noi)=res_vx;
            obj.residual_history_vy(obj.sim.noi)=res_vr;
            obj.residual_history_p(obj.sim.noi)=res_p;

            % Checking convergence
            converged = obj.has_converged(res_vx,obj.tol_vel) && obj.has_converged(res_vr,obj.tol_vel) && obj.has_converged(res_p,obj.tol_pressure);
        end
        
        function res = calculate_residual(obj, coeff, field_name)
            % For the velocity in the X and the Y direction residuals of half-control-volume boundary cells aren't taken into account
            switch field_name
                case "vx"
                    % A column and a row of zeros used to pad the data in order to make sure that there is no index out of bounds error
                    col0 = zeros(obj.sz(1)-1,1);
                    row0 = zeros(1,obj.sz(2));

                    % r = aW*fW+aE*fE+aS*fS+aN*fN+b-aP*fP
                    res_vector = [row0; coeff(2:end-1,:,1).*obj.vx(1:end-2,:); row0] + [row0; coeff(2:end-1,:,2).*obj.vx(3:end,:); row0] ... 
                        + [row0; [col0 coeff(2:end-1,2:end,3).*obj.vx(2:end-1,1:end-1)]; row0] + [row0; [coeff(2:end-1,1:end-1,4).*obj.vx(2:end-1,2:end) col0]; row0] ...
                        + [row0; coeff(2:end-1,:,6); row0] - [row0; coeff(2:end-1,:,5).*obj.vx(2:end-1,:); row0];
                    
                    % s_f = || [aP1*fP1; aP2*fP2; aP3*fP3; .... aN-1*fPN-1; aPN*fPN] ||
                    scaling_factor = norm(coeff(2:end-1,:,5).*obj.vx(2:end-1,:),1);
                    
                    % res = || r || / s_f
                    res = norm(res_vector,1)/scaling_factor;
                case "vr"
                    % A column and a row of zeros used to pad the data in order to make sure that there is no index out of bounds error
                    col0 = zeros(obj.sz(1),1);
                    row0 = zeros(1,obj.sz(2)-1);

                    % r = aW*fW+aE*fE+aS*fS+aN*fN+b-aP*fP
                    res_vector = [col0 [row0; coeff(2:end,2:end-1,1).*obj.vr(1:end-1,2:end-1)] col0] + [col0 [coeff(1:end-1,2:end-1,2).*obj.vr(2:end,2:end-1); row0] col0] ... 
                        + [col0 coeff(:,2:end-1,3).*obj.vr(:,1:end-2) col0] + [col0 coeff(:,2:end-1,4).*obj.vr(:,3:end) col0] ...
                        + [col0 coeff(:,2:end-1,6) col0] - [col0 coeff(:,2:end-1,5).*obj.vr(:,2:end-1) col0];
                    
                    % s_f = || [aP1*fP1; aP2*fP2; aP3*fP3; .... aN-1*fPN-1; aPN*fPN] ||
                    scaling_factor = norm(coeff(:,2:end-1,5).*obj.vr(:,2:end-1),1);
                    
                    % res = || r || / s_f
                    res = norm(res_vector,1)/scaling_factor;
                case "p"
                    % As it the pressure correction equation it is expected 
                    % that as the solution converges the zero field will be
                    % a solution to the pressure correction equation

                    % Hence the only term below is the b term
                    res_vector = coeff(:,:,6);
                    % s_f = || [aP1*fP1; aP2*fP2; aP3*fP3; .... aPN-1*fPN-1; aPN*fPN] ||
                    scaling_factor = norm(coeff(:,:,5).*obj.p,1);
                    
                    % res = || r || / s_f
                    res = norm(res_vector,1)/scaling_factor;    
            end
        end
        function converged = has_converged(obj, res, tol)
            converged = res<tol;
        end
        function delete_residual_history(obj)
            obj.residual_history = zeros(10000,1);
        end
    end
end