function [x_ISO_rot_state,x_ISO_ine_state,et_intersection,t_ISO,idx_inertial] = ...
    get_ISO_ine_and_rot(nodal_point_Rotating,nodal_point_Inertial, ...
    idx_ISO,options)

idx_nodal_point = nodal_point_Rotating(idx_ISO).idx;
idx_inertial    = find([nodal_point_Inertial.idx] == idx_nodal_point);
idx_inertial    = idx_inertial(1);
state_nodal_point_inertial = nodal_point_Inertial(idx_inertial).state;
et_intersection = nodal_point_Inertial(idx_inertial).et_intersection;

[t_ISO,x_ISO_ine] = ode45(@TwoBP_der,-[0 nodal_point_Inertial(idx_inertial).tof_nd],...
    state_nodal_point_inertial,options,1);
x_ISO_ine_state = flip(x_ISO_ine',2);

R_ine2rot_ISO  = cspice_sxform('ECLIPJ2000','R3BP_SE',(et_intersection + flip(t_ISO))');
for j = 1:size(R_ine2rot_ISO,3)
    x_ISO_rot_state(:,j) = R_ine2rot_ISO(:,:,j)*x_ISO_ine_state(:,j);
    x_ISO_rot_state(4,j) = x_ISO_rot_state(4,j) + x_ISO_rot_state(2,j);
    x_ISO_rot_state(5,j) = x_ISO_rot_state(5,j) - x_ISO_rot_state(1,j);
end


end