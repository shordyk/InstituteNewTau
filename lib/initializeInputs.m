% Used to initialize the input files for main runner

%% Color Packages
load Dawn.mat;
icey = cbrewer('div','BrBG',48);

%% Load Grid
load(mapFile);

%% Load Lake Locations, make polygons in polar stereographic projection
% Sl = shaperead('tc-11-451-2017-supplement/Thw_lakes_supplemental_data/Thw_lakes_outlines.shp');
% Sn = size(Sl,1);
% for i = 1:Sn
%     [x,y] = ll2ps(Sl(i).Y(~isnan(Sl(i).Y)),Sl(i).X(~isnan(Sl(i).X)));
%     Sl(i).X = x;
%     Sl(i).Y = y;
%     clear x y
% end
% clear i

%% Get speed from measures for BCs
[spd_BC_u, spd_BC_v] = measures_interp('velocity',xy(nw_bound,1),xy(nw_bound,2));
[spd_BC_u2, spd_BC_v2] = measures_interp('velocity',xy(se_bound,1),xy(se_bound,2));
[spd_BC_uL, spd_BC_vL] = measures_interp('velocity',xy(sw_bound,1),xy(sw_bound,2));
%% Load accumlation and surface temp data
[Acc, T_s] = loadALBMAP();
T = T_s(xy(:,1),xy(:,2));
