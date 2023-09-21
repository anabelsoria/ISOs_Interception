function load_lib_paths

% Add folders to path
src_Path = 'src\';
addpath(genpath(src_Path))

data_Path = 'data\';
addpath(genpath(data_Path))

dllDirectory_Path = 'src\Lambert_Russell\matlabInterface\lib\'; 
iflag = ivLam_initializeDLL(dllDirectory_Path);

% Load 3BP Rotating Frames
cspice_furnsh('src\SPICE_files\r3bp_earth-moon.tf')
cspice_furnsh('src\SPICE_files\r3bp_sun-earth.tf')
cspice_furnsh('src\SPICE_files\rssd0002.tf')
% Load Solar System Ephemeris and parameters
cspice_furnsh('src\SPICE_files\de440s.bsp')
cspice_furnsh('src\SPICE_files\naif0012.tls')
% Eros Comet ephemeris
cspice_furnsh('src\SPICE_files\2000433.bsp')
% L2 ephemeris
cspice_furnsh('src\SPICE_files\L2_de431.bsp')


end