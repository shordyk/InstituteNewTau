% Trouble shooting the model results and observational speed


% Observational
dx2 = 2e3;
smth = 4e3;
oxi = xmin-dx2*overgrab:dx2:xmax+dx2*overgrab;
oyi = ymin-dx2*overgrab:dx2:ymax+dx2*overgrab;
[oXi,oYi] = meshgrid(oxi,oyi);


b_raw =  bedmachine_interp('bed',oXi,oYi);
sf_raw =  bedmachine_interp('surface',oXi,oYi);
phi_raw = rho/rho_w.*sf_raw + (rho_w-rho)/rho_w.*b_raw;
[uO, vO] = measures_interp('velocity',oXi,oYi);
spd2  = measures_interp('speed',oXi,oYi);
overgrab = 0;
xmax = -7e5;
xmin = -12e5;
ymax =  5e5;
ymin =  -1e5;
ss2   = zeros(size(uO));

[oux ,  ouy] = gradient(uO,dx2,dx2);
[ovx ,  ovy] = gradient(vO,dx2,dx2);

oe_eff = sqrt(.5*(oux.^2 + ovy.^2) + (.5*(ouy + ovx)).^2);
[oe_effx, oe_effy] = gradient(oe_eff.^(1/3-1),dx2,dx2);

% Model speeds
% Get u,v and xy
load data_strainMesh035ISSM_centerPaulsBase.mat
dx = 259.6;
dy = 252.23 ;
xi = min(xy(:,1))-dx:dx:max(xy(:,1))+dx;
yi = (min(xy(:,2))-dy:dy:max(xy(:,2))+dy)';
%xi = xy(:,1);
%yi = xy(:,2);
[xxx,yyy] = ndgrid(xi,yi);
[Xi,Yi] = meshgrid(xi,yi);

us = scatteredInterpolant(xy(:,1),xy(:,2),u);
vs = scatteredInterpolant(xy(:,1),xy(:,2),v);


vv = vs(xxx,yyy)*3.154E7;
uu = us(xxx,yyy)*3.154E7;
ss   = zeros(size(uu));
spd = sqrt(uu.^2 + vv.^2);;

[ux ,  uy] = gradient(uu,dx,dx);
[vx ,  vy] = gradient(vv,dx,dx);

e_eff = sqrt(.5*(ux.^2 + vy.^2) + (.5*(uy + vx)).^2);
[e_effx, e_effy] = gradient(e_eff.^(1/3-1),dx,dx);

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
