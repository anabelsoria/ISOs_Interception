
%% Get the UserData
% ------------------
NoofFam = 1;
NoofEqPoints = 2;

for PointLoc = 1:NoofEqPoints
system = "sun-earth";
family = "halo"; % halo,vertical,axial,lyapunov,longp,short,butterfly
                 % dragonfly,resonant,dro,dpo,lpo
libr   = num2str(PointLoc);
branch = 'S';
T_desired_d = 140; % Days
T_desired = T_desired_d *24*60*60; % Seconds

% Use JPL Initial Conditions to propagate the orbit
url = createURL(system, family, libr,branch);
JPL = webread(url);
data = JPL.data;

% TU LU MU
TU = JPL.system.tunit;
LU = JPL.system.lunit;
MU = str2double(JPL.system.mass_ratio);

Lagrange.L1 = JPL.system.L1;
Lagrange.L2 = JPL.system.L2;
Lagrange.L3 = JPL.system.L3;
Lagrange.L4 = JPL.system.L4;
Lagrange.L5 = JPL.system.L5;

G_var = GlobalData(MU, LU, TU,Lagrange);

%% Get all the orbit parameters
% -------------------------------------------------------------------------
% Find IC of orbit with closest period
for i = 1:length(data)
    diff_per(i) = abs(cell2mat(cellfun(@str2num,data{i,1}(8),'un',0)) - T_desired/TU);
end
[~,idx] = min(diff_per);

X0(:,1,PointLoc)    = cell2mat(cellfun(@str2num,data{idx,1}(1:6),'un',0));
tPeriod(1,PointLoc) = str2num(cell2mat(data{idx,1}(8))); 

end 

HaloOrb = LyapOrbitParameters(X0,NoofFam,NoofEqPoints,tPeriod,G_var);

pathname = fullfile(cd);
%use that when you save
matfile  = fullfile(pathname, ['G_var_JPL_SE_' num2str(T_desired_d) '.mat']);
lyapfile = fullfile(pathname, ['LyapOrb_JPL_SE_' num2str(T_desired_d) '_Halo.mat']);

save(matfile,'G_var')
save(lyapfile,'HaloOrb')
