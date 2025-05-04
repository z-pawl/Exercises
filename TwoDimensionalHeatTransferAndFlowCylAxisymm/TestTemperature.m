sz=[50 500];

dx=ones(sz(1),1)/sz(1);
dr=ones(1,sz(2))/sz(2);


grid = Grid2D(dx, dr, 1e-10);
velocity = Velocity2DSIMPLE(sz,zeros(sz+[1 0]), zeros(sz+[0 1]), zeros(sz), ...
    @(x) ones(sz), @(x) ones(sz)*0, ...
    @(x) ones(sz+[1 0])*0, @(x) ones(sz+[1 0])*0, ...
    @(x) ones(sz+[0 1])*0, @(x) ones(sz+[0 1])*0, ...
    1e-7, 1e-7, 0.9, 0.3);

temp = TemperatureCondConv2DSteady(sz,230*ones(sz),@(x) ones(sz), @(x) ones(sz), @(x) ones(sz)*500, ...
    @(x) ones(sz)*0, 1e-10, 1);

temp.grid=grid;
temp.velocity=velocity;

for i=1:sz(1)
    temp.set_boundary_conditions("south",i,[2,0,0]);
    temp.set_boundary_conditions("north",i,[1,300,0]);
end
for i=1:sz(2)
    temp.set_boundary_conditions("west",i,[2,0,0]);
    temp.set_boundary_conditions("east",i,[2,0,0]);
end

sim=Simulation2D;
sim.noi=0;
temp.sim=sim;

a=false;
while ~a
    sim.noi=sim.noi+1;
    a=temp.iterate();
    if mod(sim.noi,3)==0
        fprintf('Current iteration: %d \n',sim.noi);
        fprintf('The residual is equal to: %d \n',temp.residual_history(sim.noi));
    end
end