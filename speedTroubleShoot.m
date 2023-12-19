% Trouble shooting the model results and observational speed
% Currently, want to trouble shoot grid formation. Make sure speeds match
% and tau looks right 

%% Observational

overgrab = 0;
xmax = -7e5;
xmin = -12e5;
ymax =  5e5;
ymin =  -1e5;


dx2 = 2e3;
smth = 4e3;
oxi = xmin-dx2*overgrab:dx2:xmax+dx2*overgrab;
oyi = ymin-dx2*overgrab:dx2:ymax+dx2*overgrab;
[oXi,oYi] = meshgrid(oxi,oyi);

b_raw =  bedmachine_interp('bed',oXi,oYi);
sf_raw =  bedmachine_interp('surface',oXi,oYi);
%phi_raw = rho/rho_w.*sf_raw + (rho_w-rho)/rho_w.*b_raw;
[uO, vO] = measures_interp('velocity',oXi,oYi);
spd2  = measures_interp('speed',oXi,oYi);
ss2   = zeros(size(uO));

[oux ,  ouy] = gradient(uO,dx2,dx2);
[ovx ,  ovy] = gradient(vO,dx2,dx2);

oe_eff = sqrt(.5*(oux.^2 + ovy.^2) + (.5*(ouy + ovx)).^2);
%[oe_effx, oe_effy] = gradient(oe_eff.^(1/3-1),dx2,dx2);
oe_eff_powered = oe_eff.^(1/3-1);
[oe_effx,oe_effy] = gradient(oe_eff_powered,dx2,dx2);

%% Model speeds
% Get u,v and xy
load data_strainMesh035ISSM_centerPaulsBase.mat
dx = 100;
dy = 100 ;
xi = min(xy(:,1))-dx:dx:max(xy(:,1))+dx;
yi = (min(xy(:,2))-dy:dy:max(xy(:,2))+dy);
%xi = xy(:,1);
%yi = xy(:,2);
[xxx,yyy] = ndgrid(xi,yi);
[Xi,Yi] = meshgrid(xi,yi);

us = scatteredInterpolant(xy(:,1),xy(:,2),u, 'linear', 'none');
vs = scatteredInterpolant(xy(:,1),xy(:,2),v, 'linear', 'none');


vv = vs(xxx,yyy)*3.154E7;
uu = us(xxx,yyy)*3.154E7;
ss   = zeros(size(uu));
spd = sqrt(uu.^2 + vv.^2);

[ux ,  uy] = gradient(uu,dx,dx);
[vx ,  vy] = gradient(vv,dx,dx);

e_eff = sqrt(.5*(ux.^2 + vy.^2) + (.5*(uy + vx)).^2);
%Finite Diff 
e_effx_manual = zeros(size(e_eff));
e_effy_manual = zeros(size(e_eff));

for i = 2:size(e_eff, 2) - 1
    for j = 2:size(e_eff, 1) - 1
        e_effx_manual(j, i) = (e_eff(j, i + 1) - e_eff(j, i - 1)) / (2 * dx);
        e_effy_manual(j, i) = (e_eff(j + 1, i) - e_eff(j - 1, i)) / (2 * dy);
    end
end


%e_eff_powered = e_eff.^(1/3-1);
%[e_effx, e_effy] = gradient(e_eff_powered,dx,dy);

figure
p = surf(xxx,yyy,zeros(size(ss)),spd);
hold on 
% quiver(Xi,Yi,utmp,vtmp)
title('Ice Surface Speed - Model')
set(p, 'edgecolor', 'none');
view(2)
axis equal
setFontSize(16);
c = colorbar;
c.Label.String = 'Speed [m/yr]';


figure
p = surf(oXi,oYi,zeros(size(ss2)),spd2);
hold on 
% quiver(Xi,Yi,utmp,vtmp)
title('Ice Surface Speed')
set(p, 'edgecolor', 'none');
view(2)
axis equal
setFontSize(16);
c = colorbar;
c.Label.String = 'Speed [m/yr]';

figure
p = surf(xxx,yyy,zeros(size(ss)),e_eff);
hold on 
title('model e_eff')
set(p, 'edgecolor', 'none');
view(2)
axis equal
setFontSize(16);
c = colorbar;

figure
p = surf(oXi,oYi,zeros(size(ss2)),oe_eff);
hold on 
title('obsv e_eff')
set(p, 'edgecolor', 'none');
view(2)
axis equal
setFontSize(16);
c = colorbar;
%c.Label.String = 'Speed [m/yr]';

figure
p = surf(xxx,yyy,zeros(size(uu)),e_effx_manual);
hold on 
title('model e_effx')
set(p, 'edgecolor', 'none');
view(2)
axis equal
setFontSize(16);
c = colorbar;


figure
p = surf(oXi,oYi,zeros(size(ss2)),oe_effx);
hold on 
title('obsv e_effx')
set(p, 'edgecolor', 'none');
view(2)
axis equal
setFontSize(16);
c = colorbar;
%% Debugging Grid Formulation with Tau 

xii   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
yii   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
[Xii, Yii] = ndgrid(xii-307200, yii-307200); 

% xi and yi form xxx and yyy, make square grid
tau_c = defineTau("ISSM_center");
newtau = tau_c(xy(:,1),xy(:,2),u,v)./norms([u,v],2,2);
tau_interp = scatteredInterpolant(xy(:,1),xy(:,2),newtau);

% from visualize tau
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




figure
p = surf(xxx,yyy,zeros(size(ss)),tau_interp(xxx, yyy), "EdgeColor", "none");
hold on
%contour(xi,yi,spd, [30, 30] , 'k--','HandleVisibility','off')
%contour(xi,yi,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
%contour(xi,yi,spd, [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','color',rgb('gray'),'linewidth',2)
title('Tau - ISSM Center')
colorbar
view(2)

[XX,YY] = ndgrid(xy(:,1),xy(:,2));

figure
p = surf(xxx,yyy,zeros(size(ss)),tau_interp(XX, YY), "EdgeColor", "none");
hold on
%contour(xi,yi,spd, [30, 30] , 'k--','HandleVisibility','off')
%contour(xi,yi,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
%contour(xi,yi,spd, [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','color',rgb('gray'),'linewidth',2)
title('Tau - ISSM Center XX')
colorbar
view(2)





