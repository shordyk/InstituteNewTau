% Used to initialize the input files for main runner

%% Confirm that ISMIP data is downloaded, otherwise downloads and unpacks
if(exist("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc",'file')~=2)
    warning('Missing initMIP data, downloading now');
    system('bash lib/downloadIsMipData.sh');
    disp('Download Finished') ;
end

%% Color Packages
load Dawn.mat;
icey = cbrewer('div','BrBG',48);

%% Load Grid
load(mapFile);

%% Get speed from measures for BCs
[spd_BC_u, spd_BC_v] = measures_interp('velocity',xy(nw_bound,1),xy(nw_bound,2));
[spd_BC_u2, spd_BC_v2] = measures_interp('velocity',xy(se_bound,1),xy(se_bound,2));
[spd_BC_uL, spd_BC_vL] = measures_interp('velocity',xy(sw_bound,1),xy(sw_bound,2));

%% Load accumlation and surface temp data
% These are large data sets and are scrubbed before saving to data file
[Acc, T_s] = loadALBMAP();
T = T_s(xy(:,1),xy(:,2));
