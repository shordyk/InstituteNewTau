% Open and work with ALPMAP data for SMB and Temp data
clear
Tm  = 273;     %Melt point [K]
rho = 917;     %Ice density[kg/m^3]
Cp  = 2050;    %specific heat of ice [J/Kg/K]
K   = 2.1;     %thermal conductivity of ice [W/m/K]
A   = 2.4e-24; %prefactor [Pa^-3 s^-1]
n   = 3;       %Glens law power
%% Load Data
struct = ncinfo('ALBMAPv1.nc');
xi = ncread('ALBMAPv1.nc','x1');
yi = ncread('ALBMAPv1.nc','y1');
temp = ncread('ALBMAPv1.nc','temp'); 
acca = ncread('ALBMAPv1.nc','acca'); 
accr = ncread('ALBMAPv1.nc','accr'); 

%% Load Grid
load('../workingGrid4.mat')
[Xi,Yi] = ndgrid(xi,yi);
triTemp = griddedInterpolant(Xi,Yi,temp);
triACCA = griddedInterpolant(Xi,Yi,acca);
triACCR = griddedInterpolant(Xi,Yi,accr);
 

%% Plot
figure
    trisurf(t,xy(:,1),xy(:,2),triTemp(xy(:,1),xy(:,2)),'edgecolor','none');
    colorbar
    view(2)
    title('average surface temperature')

figure
    subplot(211)
        trisurf(t,xy(:,1),xy(:,2),triACCA(xy(:,1),xy(:,2)),'edgecolor','none');
        colorbar
        view(2)
        title('ACCA')
    subplot(212)
        trisurf(t,xy(:,1),xy(:,2),triACCR(xy(:,1),xy(:,2)),'edgecolor','none');
        colorbar
        view(2)
        title('ACCR')
