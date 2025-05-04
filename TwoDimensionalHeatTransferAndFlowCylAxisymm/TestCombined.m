sz=[450 75];

dx=0.3*ones(sz(1),1)/sz(1);
dr=0.05*ones(1,sz(2))/sz(2);


grid = Grid2D(dx, dr, 1e-10);


vel=Velocity2DSIMPLE(sz,0.0*ones(sz+[1 0]),zeros(sz+[0 1]),zeros(sz),@(x) 1000*ones(sz), @(x) 1e-3*ones(sz), ...
    @(x) zeros(sz+[1 0]), @(x) zeros(sz+[1 0]), @(x) zeros(sz+[0 1]), @(x) zeros(sz+[0 1]), 1e-7, 1e-7, 0.85, 0.3);

vel.grid=grid;
sim=Simulation2D;
sim.noi=0;
vel.sim=sim;

for i=1:sz(1)
    vel.set_boundary_conditions("vr","south",i,[2,0,0]);
    vel.set_boundary_conditions("vr","north",i,[1,0,0]);

    vel.set_boundary_conditions("p","south",i,[2,0,0]);
    vel.set_boundary_conditions("p","north",i,[2,0,0]);
end
for i=1:sz(2)
    vel.set_boundary_conditions("vx","west",i,[1,0.01,0]);
    vel.set_boundary_conditions("vx","east",i,[2,0,0]);

    vel.set_boundary_conditions("p","west",i,[2,0,0]);
    vel.set_boundary_conditions("p","east",i,[1,0,0]);
end
% vel.set_boundary_conditions("p","west",1,[1,0,0]);
for i=1:sz(1)+1
    vel.set_boundary_conditions("vx","south",i,[1,0,0]);
    vel.set_boundary_conditions("vx","north",i,[1,0,0]);
end
for i=1:sz(2)+1
    vel.set_boundary_conditions("vr","west",i,[1,0,0]);
    vel.set_boundary_conditions("vr","east",i,[2,0,0]);
end


temp = TemperatureCondConv2DSteady(sz,800*ones(sz),@(x) 47.747*ones(sz), @(x) 50*ones(sz), ...
    @(x) zeros(sz), @(x) zeros(sz), 1e-7, 1);

temp.grid=grid;
temp.velocity=vel;
temp.sim=sim;

for i=1:sz(1)
    temp.set_boundary_conditions("north",i,[1,900,0]);
    temp.set_boundary_conditions("south",i,[2,0,0]);
end
for i=1:sz(2)
    temp.set_boundary_conditions("west",i,[1,700,0]);
    temp.set_boundary_conditions("east",i,[2,0,0]);
end

vel.temp=temp;
%vel.srcy_f=@(x) boussinesq_approximation_y(x);

a=false;
while ~a
    sim.noi=sim.noi+1;
    v=vel.iterate();
    t=temp.iterate();
    a=v && t;
    if mod(sim.noi,10)==0
        fprintf('Current iteration: %d \n',sim.noi);
        fprintf('The vx residual is equal to: %d \n',vel.residual_history_vx(sim.noi));
        fprintf('The vr residual is equal to: %d \n',vel.residual_history_vy(sim.noi));
        fprintf('The p residual is equal to: %d \n',vel.residual_history_p(sim.noi));
        fprintf('The T residual is equal to: %d \n',temp.residual_history(sim.noi));
    end

    if sim.noi==10000
        break
    end
end



% Plotting the results
% Grid
[X,Y]=meshgrid(grid.pos_cent_x,grid.pos_cent_r);

% Interpolated velocities
vx_int=(vel.vx(1:end-1,:)+vel.vx(2:end,:))/2;
vr_int=(vel.vr(:,1:end-1)+vel.vr(:,2:end))/2;
v_magn_int=sqrt(vx_int.^2+vr_int.^2);


% Values have to be transposed - probably because MATLAB uses column major order
% Temperature graph
figure(1);
colormap(jet);
contourf(X,Y,transpose(temp.temp),30);

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
semilogy(vel.residual_history_vx(1:sim.noi),'DisplayName','vx residual');
hold on;
semilogy(vel.residual_history_vy(1:sim.noi),'DisplayName','vr residual');
semilogy(vel.residual_history_p(1:sim.noi),'DisplayName','p residual');
semilogy(temp.residual_history(1:sim.noi),'DisplayName','T residual');
