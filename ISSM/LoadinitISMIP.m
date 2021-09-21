clear;
% https://tc.copernicus.org/articles/13/1441/2019/tc-13-1441-2019.pdf
%%
% filename = "JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc";
% filename_init = "JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc";

% filename = "ARC_PISM1_ctrl/strbasemag_AIS_ARC_PISM1_ctrl.nc";
filename2 = "IGE_ELMER_ctrl/strbasemag_AIS_IGE_ELMER_ctrl.nc";

% 
filename = "JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc";
x   = ncread("JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
y   = ncread("JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
tau = ncread("JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
% x   = ncread(filename,"x");
% y   = ncread(filename,"y");
% tau = ncread(filename,"strbasemag");

tau2 = ncread(filename2,"strbasemag");
x2   = ncread(filename2,"x");
y2   = ncread(filename2,"y");

%% Clean some values
x = x - 3072000;
y = y - 3072000;
x2 = x2;
y2 = y2;
tau(isnan(tau)) = 0;
tau2(isnan(tau2)) = 0;

%%
[xx,yy] = ndgrid(x,y);
[xx2,yy2] = ndgrid(x2,y2);
uB = griddedInterpolant(xx,yy,tau(:,:,21));
uB2 = griddedInterpolant(xx2,yy2,tau2(:,:,21));

% figure(1)
% clf
% surf(x,y,tau'/1e3,'edgecolor','none')
% caxis([0 150]);
% colorbar
% view(2)
% 
x0 = mean(x);
y0 = mean(y);

figure(1)
clf
surf(x,y,uB(xx,yy)'/1e3,'edgecolor','none')
hold on
scatter3(x0,y0,400,[],'k','filled');
caxis([0 150]);
colorbar
view(2)
tmp1 = split(filename,"/");
title(tmp1(1))

%%
figure(2)
clf
surf(x2,y2,uB2(xx,yy)'/1e3,'edgecolor','none')
hold on
caxis([0 150]);
colorbar
view(2)
tmp2 = split(filename2,"/");
title(tmp2(1))
%%
figure(3)
clf
surf(x2,y2,(uB(xx,yy) - uB2(xx,yy))'/1e3,'edgecolor','none')
hold on
caxis([-150 150]);
colorbar
colormap redblue
title(tmp1(1) + " - " + tmp2(1))
view(2)
%%
load ../gridInstitute24000.mat
figure(4)
clf
subplot(131)
trisurf(t,xy(:,1),xy(:,2),uB(xy(:,1),xy(:,2))/1e3,...
       'edgecolor','none')
title(tmp1(1))
view(2)
colorbar

subplot(132)
trisurf(t,xy(:,1),xy(:,2),uB2(xy(:,1),xy(:,2))/1e3,...
       'edgecolor','none')
title(tmp2(1))
view(2)
colorbar

ax(3) = subplot(133);
trisurf(t,xy(:,1),xy(:,2),(uB(xy(:,1),xy(:,2)) - uB2(xy(:,1),xy(:,2)))/1e3,...
       'edgecolor','none')
colormap(ax(3), redblue)
title(tmp1(1) + " - " + tmp2(1))
view(2)
colorbar