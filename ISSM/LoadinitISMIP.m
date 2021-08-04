clear;
% https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.pdf

x   = ncread("JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
y   = ncread("JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
tau = ncread("JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","strbasemag");
tau2 = ncread("JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
x2   = ncread("JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","x");


x = x - 3072000;
y = y - 3072000;

[xx,yy] = ndgrid(x,y);
uB = griddedInterpolant(xx,yy,tau');
uB2 = griddedInterpolant(xx,yy,tau2(:,:,21)');

figure(1)
clf
surf(x,y,tau'/1e3,'edgecolor','none')
caxis([0 150]);
colorbar
view(2)

x0 = mean(x);
y0 = mean(y);

figure(2)
clf
surf(x,y,uB(xx,yy)/1e3,'edgecolor','none')
hold on
scatter3(x0,y0,400,[],'k','filled');
caxis([0 150]);
colorbar
view(2)

%%
figure(3)
clf
surf(x,y,(uB(xx,yy) - uB2(xx,yy))/1e3,'edgecolor','none')
hold on
scatter3(x0,y0,400,[],'k','filled');
% caxis([0 150]);
colorbar
view(2)