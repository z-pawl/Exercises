sz=[75 75];

dx=ones(sz(1),1)/sz(1);
dy=ones(1,sz(2))/sz(2);


grid = Grid2D(dx, dy);

Re=1000;

function buyoancy_term_y = boussinesq_approximation_y(vel)
    rho0=1;
    g=-9.81;
    ref_temp=800;
    sz = [75 75];
    col0 = zeros(sz(1),1);
    
    temp=vel.temp.temp;
    temp_int=([col0 temp].*[vel.grid.dy 1e25]+[temp col0].*[1e25 vel.grid.dy])./([vel.grid.dy 1e25]+[1e25 vel.grid.dy]);

    %buyoancy_term_y=zeros(sz+[0 1]);
    buyoancy_term_y=(2*ones(sz+[0 1])-temp_int/ref_temp)*rho0*g;

end

vel=Velocity2DSIMPLE(sz,0.0*ones(sz+[1 0]),zeros(sz+[0 1]),zeros(sz),@(x) ones(sz), @(x) 1/Re*ones(sz), ...
    @(x) zeros(sz+[1 0]), @(x) zeros(sz+[1 0]), @(x) zeros(sz+[0 1]), @(x) zeros(sz+[0 1]), 1e-7, 1e-7, 0.55, 0.2);

vel.grid=grid;
sim=Simulation2D;
sim.noi=0;
vel.sim=sim;

for i=1:sz(1)
    vel.set_boundary_conditions("vy","south",i,[1,0,0]);
    vel.set_boundary_conditions("vy","north",i,[1,0,0]);

    vel.set_boundary_conditions("p","south",i,[2,0,0]);
    vel.set_boundary_conditions("p","north",i,[2,0,0]);
end
for i=1:sz(2)
    vel.set_boundary_conditions("vx","west",i,[1,0,0]);
    vel.set_boundary_conditions("vx","east",i,[1,0,0]);

    vel.set_boundary_conditions("p","west",i,[2,0,0]);
    vel.set_boundary_conditions("p","east",i,[2,0,0]);
end
vel.set_boundary_conditions("p","west",1,[1,0,0]);
for i=1:sz(1)+1
    vel.set_boundary_conditions("vx","south",i,[1,0,0]);
    vel.set_boundary_conditions("vx","north",i,[1,0,0]);
end
for i=1:sz(2)+1
    vel.set_boundary_conditions("vy","west",i,[1,0,0]);
    vel.set_boundary_conditions("vy","east",i,[1,0,0]);
end


temp = TemperatureCondConv2DSteady(sz,800*ones(sz),@(x) 407.747*ones(sz), @(x) ones(sz), ...
    @(x) zeros(sz), @(x) zeros(sz), 1e-7, 1);

temp.grid=grid;
temp.velocity=vel;
temp.sim=sim;

for i=1:sz(1)
    temp.set_boundary_conditions("north",i,[2,0,0]);
    temp.set_boundary_conditions("south",i,[2,0,0]);
end
for i=1:sz(2)
    temp.set_boundary_conditions("west",i,[1,700,0]);
    temp.set_boundary_conditions("east",i,[1,900,0]);
end

vel.temp=temp;
vel.srcy_f=@(x) boussinesq_approximation_y(x);

a=false;
while ~a
    sim.noi=sim.noi+1;
    v=vel.iterate();
    t=temp.iterate();
    a=v && t;
    if mod(sim.noi,10)==0
        fprintf('Current iteration: %d \n',sim.noi);
        fprintf('The vx residual is equal to: %d \n',vel.residual_history_vx(sim.noi));
        fprintf('The vy residual is equal to: %d \n',vel.residual_history_vy(sim.noi));
        fprintf('The p residual is equal to: %d \n',vel.residual_history_p(sim.noi));
        fprintf('The T residual is equal to: %d \n',temp.residual_history(sim.noi));
    end

    if sim.noi==10000
        break
    end
end



% Plotting the results
% Grid
[X,Y]=meshgrid(grid.pos_cent_x,grid.pos_cent_y);

% Interpolated velocities
vx_int=(vel.vx(1:end-1,:)+vel.vx(2:end,:))/2;
vy_int=(vel.vy(:,1:end-1)+vel.vy(:,2:end))/2;
v_magn_int=sqrt(vx_int.^2+vy_int.^2);


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
contourf(X,Y,transpose(vy_int),30);

% Streamlines
figure(5);
streamslice(X,Y,transpose(vx_int),transpose(vy_int),2);

% Residuals
figure(6);
semilogy(vel.residual_history_vx(1:sim.noi),'DisplayName','vx residual');
hold on;
semilogy(vel.residual_history_vy(1:sim.noi),'DisplayName','vy residual');
semilogy(vel.residual_history_p(1:sim.noi),'DisplayName','p residual');
semilogy(temp.residual_history(1:sim.noi),'DisplayName','T residual');

