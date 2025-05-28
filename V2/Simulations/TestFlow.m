offset_R = 1e-5;
delta_R = 1;
L = 1;

sz = [5 5];

dx = L * ones(sz(1),1) / sz(1);
dr = delta_R * ones(1,sz(2)) / sz(2);

grid = Grid2D(dx, dr, offset_R);

flow = FlowComponent(grid, zeros(sz), zeros(sz), zeros(sz), 1, @(x) ones(sz), @(x) ones(sz), @(x) zeros(sz), @(x) zeros(sz), @(x) zeros(sz), @(x) zeros(sz), 1e-6, 1e-6, 0.55, 0.2, 1, 1, 30);

% Boundary conditions
% West
for j = 1:sz(2)
    flow.vx_bds.set_boundary_condition("x", [1 j], 1, 1, 0, 0.2);
    flow.vr_bds.set_boundary_condition("x", [1 j], 1, 1, 0, 0);
    flow.p_bds.set_boundary_condition("x", [1 j], 1, 0, 1, 0);
end
% East
for j = 1:sz(2)
    flow.vx_bds.set_boundary_condition("x", [sz(1)+1 j], -1, 0, 1, 0);
    flow.vr_bds.set_boundary_condition("x", [sz(1)+1 j], -1, 0, 1, 0);
    flow.p_bds.set_boundary_condition("x", [sz(1)+1 j], -1, 1, 0, 0);
end
% South
for i = 1:sz(1)
    flow.vx_bds.set_boundary_condition("r", [i 1], 1, 0, 1, 0);
    flow.vr_bds.set_boundary_condition("r", [i 1], 1, 0, 1, 0);
    flow.p_bds.set_boundary_condition("r", [i 1], 1, 0, 1, 0);
end
% North
for i = 1:sz(1)
    flow.vx_bds.set_boundary_condition("r", [i sz(2)+1], -1, 1, 0, 0);
    flow.vr_bds.set_boundary_condition("r", [i sz(2)+1], -1, 1, 0, 0);
    flow.p_bds.set_boundary_condition("r", [i sz(2)+1], -1, 0, 1, 0);
end

flow.convert_pressure_bd_conditions();

flow.update_properties();
[coeff_vx, coeff_vr] = flow.get_coefficients_v();
flow.update_face_velocities(flow.vx, flow.vr, coeff_vx(:,:,1), coeff_vr(:,:,1));

a = false;
while ~a
    a = flow.iterate();
    if mod(flow.noi,10)==0
        fprintf('Current iteration: %d \n', flow.noi);
        fprintf('The vx residual is equal to: %d \n', flow.residual_history_vx(flow.noi));
        fprintf('The vr residual is equal to: %d \n', flow.residual_history_vr(flow.noi));
        fprintf('The p residual is equal to: %d \n', flow.residual_history_p(flow.noi));
    end
end