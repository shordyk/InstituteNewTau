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


v = vs(Xi,Yi);
u = us(Xi,Yi);
ss   = zeros(size(u));
spd = sqrt(u.^2 + v.^2);

figure
p = surf(Xi,Yi,zeros(size(ss)),spd);
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
