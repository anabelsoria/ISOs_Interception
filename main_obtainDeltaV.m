% Author:        Anabel Soria

clc
clear all
close all

load_lib_paths


%% Inputs
idx_ISO = 8;
libPtNo = 2;
orbitNo = 1;
dimension_po = 300;

%% ---------------------------LYAPUNOV-------------------------------------

Gvar_name = ['G_var_JPL_SE_' num2str(dimension_po) '.mat' ];
G_var = importdata(Gvar_name);
LyapOrb_name = ['LyapOrb_JPL_SE_' num2str(dimension_po) '.mat' ];
LyapOrb = importdata(LyapOrb_name);

[t,~,x,~,~]      = StateTransAndX(LyapOrb(libPtNo).IC,...
                    1*LyapOrb(libPtNo).time(orbitNo),G_var.Constants.mu);
x_Lyap.rot_state = x';
x_Lyap.t         = t';

%% ---------------------------ISO------------------------------------------
ISO_data = load('iso_orbits_detectdata_v3.mat');

nodal_point_Inertial = importdata('nodal_point_Inertial_Total.mat');
nodal_point_Rotating = importdata('nodal_point_Rotating_0.7_3.mat');

options = odeset('Reltol',1e-12,'Abstol',1e-12);

for i = 1
    [rot,ine,et_intersection,t_ISO] = get_ISO_ine_and_rot(nodal_point_Rotating,...
        nodal_point_Inertial,idx_ISO,options);
    x_ISO(i).ine_state = ine;
    x_ISO(i).rot_state = rot;
end

%% ---------------------------MANIFOLDS------------------------------------
Npoints = 100;
epsilon = 100;
epsilon = epsilon/G_var.Constants.l;
dir = 1;
delta_v = 0;
type = 'unstable';
switch type 
    case 'stable'
        eigVec = LyapOrb(libPtNo).Eigens(orbitNo).S_EigVec;
    case 'unstable'
        eigVec = LyapOrb(libPtNo).Eigens(orbitNo).US_EigVec;
end
if dimension_po == 175
    eps = 100;
else
    eps = 1000;
end 

LyapX0_name = ['Lyap_X0_' num2str(dimension_po) '_eps_'...
    num2str(eps) '_exterior.mat'];
load(LyapX0_name)

%% ---------------------------TRANSFER-------------------------------------
warning('off','all')
tspan_manifold = [0.000001:0.1:4]; % N x Lyap Orbit Period
tspan_lambert  = [300:2:410]; % Days [600:2:700];      % 
q0             = [0.5,0.5];        % km/s

x_arrival =  x_ISO(1).rot_state(:,end);

fsolve_opts = optimoptions('fsolve','Algorithm','Levenberg-Marquardt','Display','off',...
    'MaxIterations',50,'FunctionTolerance',1e-12);
tic
for k = 1:Npoints
    % Restart always q0
    q0 = [0.5,0.5];
    n = 0;
    for tof_lambert_0 = tspan_lambert
        n = n + 1;
        m = 0;
        for tof_manifold_0 = tspan_manifold
            m = m + 1;
            
            % Define function
            fun =  @(x)fun_fsolveDeltaV(x,Lyap_X0(k,:),tof_manifold_0,tof_lambert_0,G_var.Constants.l,G_var.Constants.T,...
                    x_arrival,LyapOrb(libPtNo).time(orbitNo),LyapOrb(libPtNo).IC(orbitNo,:),G_var.Constants.mu);
            % Call function to solve
            [q,fval,EXITFLAG,output] = fsolve(fun,q0,fsolve_opts);

            sol(k).v(m,n,:) = q;
            sol(k).tof_manifold(m,n) = tof_manifold_0;
            sol(k).tof_lambert(m,n)  = tof_lambert_0;
            sol(k).startLyap(m,n,:)  = Lyap_X0(k,:);
            sol(k).iterations(m,n)   = output.iterations;
            sol(k).funcount(m,n)     = output.funcCount;
            
            if EXITFLAG > 0 
                q0 = q;
                sol(k).DeltaV(m,n) = norm(q);
            else
                q0 = [ 0.5,0.5];
                sol(k).DeltaV(m,n) = NaN;
            end
        end
    end
end
timer =  timer + toc;
