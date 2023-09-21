function [eq,x_manifold,x_lambert] = fun_fsolveDeltaV(q,Lyap_X0,tof_manifold,tof_lambert,LU,...
    TU,x_ISO,LyapOrb_period,LyapOrb_period_IC,MU) 
%#codegen

% Adimensionalize time of flight
tof_manifold_nd = tof_manifold*LyapOrb_period;
tof_lambert_nd  = tof_lambert*24*60*60/TU;

% ------------------------------- Manifold -------------------------------

if tof_manifold > 0.000001
    options = odeset('Reltol',1e-12,'Abstol',1e-12);
    [~,x] = ode89(@(t,x)CRes3BP_EOM(t,x,MU),[0 tof_manifold_nd],Lyap_X0,options);
    x_manifold = x';
else
    x_manifold = LyapOrb_period_IC';
end

% ---------------------------- CR3BP Transfer ----------------------------

% Forward integration from last state manifold
ic_lambert = x_manifold(:,end);
% Add Delta V
deltavx = q(1)/(LU/TU);
deltavy = q(2)/(LU/TU);
ic_lambert(4) = ic_lambert(4) + deltavx;
ic_lambert(5) = ic_lambert(5) + deltavy;

% Integration with ode45 using a custom timeout mechanism

if tof_lambert > 0.000001
    opts  = odeset('Reltol',1e-12,'Abstol',1e-12,'Events', @(t,x)eventFunction(t,x,LU,MU));
    [~,x,~,~,ie] = ode45(@(t,x)CRes3BP_EOM(t,x,MU),[0 tof_lambert_nd], ic_lambert,opts);
else
    x = x_manifold(:,end)';
    ie = [];
end


if isempty(ie)
    x_lambert = x';
else
    x_lambert = [1000;1000];
end


% F(x) to solve
eq = coder.nullcopy(zeros(2, 1));
eq(1) = x_lambert(1,end) - x_ISO(1);
eq(2) = x_lambert(2,end) - x_ISO(2);

end