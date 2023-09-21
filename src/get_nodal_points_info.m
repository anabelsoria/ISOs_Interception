function nodal_point = get_nodal_points_info(ISO_data, TU, LU, m_SUN_nd,...
    rotation, min_rad,max_rad)

p = 0;
for idx_ISO = 1:length(ISO_data.activeISOs.detect_state)
    % Choose which ISO from excel
    ic_ISO         = ISO_data.activeISOs.detect_state(idx_ISO,:);
    ic_ISO_nd(1:3) = ic_ISO(1:3)/LU;
    ic_ISO_nd(4:6) = ic_ISO(4:6)/(LU/TU);
    
    % Convert time detection JD to ET
    t_detection = ISO_data.activeISOs.detect_JD; % JD
    JD_2000  = 2451545.0;
    t_detection = (t_detection - JD_2000)*24*60*60; % ET
    
    % NUMERIC PROPAGATION ISO
    T_po     = 365*5;               % 5 years in days
    tspan    = [0 T_po]*(60*60*24); % seconds
    tspan_nd = tspan/TU;            % non-dimensional
    % Stop integration if ecliptic plane reached
    options  = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events', @myEvent);
    [t,x_ISO_inertial] = ode45(@TwoBP_der,tspan_nd,ic_ISO_nd,options,m_SUN_nd);
    x_ISO_inertial = x_ISO_inertial';
    
    % Make change to rotation:
    if rotation == true
        et_intersection = t_detection + t(end)*TU;
        R_ine2rot   = cspice_sxform('ECLIPJ2000','R3BP_SE',et_intersection);
        x_ISO_rot(1:6,:) = R_ine2rot*x_ISO_inertial(1:6,end);
        last_state = x_ISO_rot;
    else
        last_state = x_ISO_inertial(1:6,end);
    end
    
    if abs(last_state(3)) < 1e-3
        r_node_SUN = vecnorm(last_state(1:2));
        if exist('min_rad') == 1
            if r_node_SUN >= min_rad && r_node_SUN <= max_rad
                p = p + 1;
                nodal_point(p).state             = last_state;
                nodal_point(p).et_intersection   = t_detection + t(end)*TU;
                nodal_point(p).idx               = idx_ISO;
                nodal_point(p).id        = ISO_data.activeISOs.ssm_id(idx_ISO);
                nodal_point(p).tof_nd            = t(end);
                nodal_point(p).type              = 'active';
            end
        else
            if p == 1 display('Saving all nodal points'); end
            p = p + 1;
            nodal_point(p).state             = last_state;
            nodal_point(p).et_intersection   = t_detection + t(end)*TU;
            nodal_point(p).idx               = idx_ISO;
            nodal_point(p).id        = ISO_data.activeISOs.ssm_id(idx_ISO);
            nodal_point(p).tof_nd            = t(end);
            nodal_point(p).type              = 'active';
        end
    end
end

%% Repeat for inert
for idx_ISO_inert = 1:length(ISO_data.inertISOs.detect_state)
    % Choose which ISO from excel
    ic_ISO         = ISO_data.inertISOs.detect_state(idx_ISO_inert,:);
    ic_ISO_nd(1:3) = ic_ISO(1:3)/LU;
    ic_ISO_nd(4:6) = ic_ISO(4:6)/(LU/TU);
    
    % Convert time detection JD to ET
    t_detection = ISO_data.inertISOs.detect_JD; % JD
    JD_2000  = 2451545.0;
    t_detection = (t_detection - JD_2000)*24*60*60; % ET
    
    % NUMERIC PROPAGATION ISO
    T_po     = 365*5;              % days
    tspan    = [0 T_po]*(60*60*24); % seconds
    tspan_nd = tspan/TU;            % non-dimensional
    options  = odeset('RelTol',1e-12,'AbsTol',1e-12,'Events', @myEvent);
    
    [t,x_ISO_inertial] = ode45(@TwoBP_der,tspan_nd,ic_ISO_nd,options,m_SUN_nd);
    x_ISO_inertial = x_ISO_inertial';
    
    % Make change to rotation:
    if rotation == true
        et_intersection = t_detection + t(end)*TU;
        R_ine2rot   = cspice_sxform('ECLIPJ2000','R3BP_SE',et_intersection);
        x_ISO_rot(1:6,:) = R_ine2rot*x_ISO_inertial(1:6,end);
        last_state = x_ISO_rot;
    else
        last_state = x_ISO_inertial(1:6,end);
    end
    
    if abs(last_state(3)) < 1e-3
        r_node_SUN = vecnorm(last_state(1:2));
        if exist('min_rad') == 1
            if r_node_SUN >= min_rad && r_node_SUN <= max_rad
                p = p + 1;
                nodal_point(p).state             = last_state;
                nodal_point(p).et_intersection   = t_detection + t(end)*TU;
                nodal_point(p).idx               = idx_ISO_inert;
                nodal_point(p).id        = ISO_data.inertISOs.ssm_id(idx_ISO_inert);
                nodal_point(p).tof_nd            = t(end);
                nodal_point(p).type              = 'inert';
            end
        else
            p = p + 1;
            nodal_point(p).state             = last_state;
            nodal_point(p).et_intersection   = t_detection + t(end)*TU;
            nodal_point(p).idx               = idx_ISO_inert;
            nodal_point(p).id        = ISO_data.inertISOs.ssm_id(idx_ISO_inert);
            nodal_point(p).tof_nd            = t(end);
            nodal_point(p).type              = 'inert';
        end
    end
end
