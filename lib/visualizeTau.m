clear
close all

warning('off','MATLAB:imagesci:netcdf:fillValueTypeMismatch'); % try not to warn here
xi   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
yi   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
tau  = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
warning('on','MATLAB:imagesci:netcdf:fillValueTypeMismatch'); 
[xx,yy] = ndgrid(xi - 3072000,yi - 3072000);

tau_interp = griddedInterpolant(xx,yy,tau(:,:,21));

%%
figure(1)
imagesc(xi,yi,tau(:,:,21)')
set(gca,'YDir','normal')
colorbar

%%
xrange = [-1000000 -825000];
yrange = [ 170000   340000];

x_small = linspace(xrange(1),xrange(2),1000);
y_small = linspace(yrange(1),yrange(2),1000);
[X_small,Y_small] = ndgrid(x_small,y_small);
spd = measures_interp('speed',X_small,Y_small);

figure(2)
surf(X_small,Y_small,zeros(size(X_small)),tau_interp(X_small,Y_small),'edgecolor','none')
hold on
contour(X_small,Y_small,spd,[10,10],'k:')
contour(X_small,Y_small,spd,[30,30],'k--')
contour(X_small,Y_small,spd,[100,300],'k-')
contour(X_small,Y_small,spd,[1000,1000],'k-','linewidth',2)
colorbar
title('Basal Strength')
view(2)


figure(3)
%surf(X_small,Y_small,zeros(size(X_small)),spd,20,'edgecolor','none')
contourf(X_small,Y_small,spd,'LineColor', 'none')
%set(gca,'YDir','normal')
title('Speed of Measures')
colorbar

