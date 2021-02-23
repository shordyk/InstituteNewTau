% Mapping strain and driving force for institute methods of Van Der Veen 1989
clear; close all
load Dawn.mat
icey = cbrewer('div','BrBG',48);
rho = 917;
rho_w = 1000;
g = 9.81;
B = 7.5e7; % A = 2.4e-24 Pa^(-3) s^(-1)
overgrab = 0;
xmax = -7e5;
xmin = -10e5;
ymax =  5e5;
ymin =  1e5;

dx = 2e3;
smth = 6e3;
xi = xmin-dx*overgrab:dx:xmax+dx*overgrab;
yi = ymin-dx*overgrab:dx:ymax+dx*overgrab;
[Xi,Yi] = meshgrid(xi,yi);

% Masking
% msk = scatteredInterpolant(xy(:,1),xy(:,2),ones(size(xy(:,1))),'nearest','none');
% mask = msk(Xi,Yi);
% mask(isnan(mask)) = 0;
 mask = ones(size(Xi)); %make mask all ones for now

%Raw fields
b_raw =  bedmachine_interp('bed',Xi,Yi);
sf_raw =  bedmachine_interp('surface',Xi,Yi);
phi_raw = rho/rho_w.*sf_raw + (rho_w-rho)/rho_w.*b_raw;

[u, v] = measures_interp('velocity',Xi,Yi);

%% Smoothing/Masking
% Gauss
b = imgaussfilt(b_raw,2);
u = imgaussfilt(u,smth/dx) / 3.154e7;
v = imgaussfilt(v,smth/dx) / 3.154e7;
spd = sqrt(u.^2 + v.^2);
sf = imgaussfilt(sf_raw,smth/dx);
phi = imgaussfilt(phi_raw,2);


% mean filtering
% hh = ones(floor(smth/dx))/(floor(smth/dx).^2);
% u = filter2(hh,u) / 3.154e7;
% v = filter2(hh,v) / 3.154e7;
% surf = filter2(hh,surf);

sf(sf < b) = b(sf < b);
spd2  = measures_interp('speed',Xi,Yi);

%Gradients
h = sf-b;
[ux ,  uy] = gradient(u,dx,dx);
[vx ,  vy] = gradient(v,dx,dx);
[sx ,  sy] = gradient(sf,dx,dx);
[spdx ,  spdy] = gradient(spd,dx,dx);
[spdxx, spdxy] = gradient(spdx,dx,dx);
[spdyx, spdyy] = gradient(spdy,dx,dx);


% Effective Strain
e_eff = sqrt(.5*(ux.^2 + vy.^2) + (.5*(uy + vx)).^2);
[e_effx, e_effy] = gradient(e_eff.^(1/3-1),dx,dx);

% Deviatoric Stress
prefactor = B * e_eff.^(1/3 - 1);

% Driving Stress
Tdx = -rho * g * h .* sx.^2;
Tdy = -rho * g * h .* sy.^2;
Td  = sqrt(Tdx.^2 +  Tdy.^2);

ss   = zeros(size(u));
angs = zeros(size(u));
dr   = zeros(size(u));
lon  = zeros(size(u));
lat  = zeros(size(u));
bed  = zeros(size(u));
h = sf-b;
for i = 2:length(xi)-1
    for j = 2:length(yi)-1
        ui = u(j,i);
        vi = v(j,i);
        ang = atan(vi/ui);
        if(ui < 0) 
            ang = ang + pi;
        end
        vv = [cos(ang), sin(ang)]; %Direction Vectors along flow
        vv_t = [-sin(ang), cos(ang)];%Direction Vectors Perp to flow
        R = [ cos(ang) sin(-ang) ; sin(ang) cos(ang) ];
        ss(j,i) = max([ui , vi] * R); %x speed in new ref frame
        dr(j,i) = -(vv(1)*sx(j,i) + vv(2)*sy(j,i))* rho * g * h(j,i); %Driving Force
        
        lon(j,i) =  2*B*((vv(1)*sx(j,i) + vv(2)*sy(j,i)) .* e_eff(j,i).^(1/3-1) .* (vv(1)*spdx(j,i) + vv(2)*spdy(j,i))...
                    + h(j,i) .* (vv(1)*e_effx(j,i) + vv(2)*e_effy(j,i)) .* (vv(1)*spdx(j,i) + vv(2)*spdy(j,i))...
                    + h(j,i) .* e_eff(j,i).^(1/3-1) .* (spdxx(j,i).*vv(1).^2 + spdxy(j,i).*vv(1).*vv(2) + spdyx(j,i).*vv(1).*vv(2) + spdyy(j,i).*vv(2).^2));
                
        lat(j,i) =  2*B*((vv_t(1)*sx(j,i) + vv_t(2)*sy(j,i)) .* e_eff(j,i).^(1/3-1) .* (vv_t(1)*spdx(j,i) + vv_t(2)*spdy(j,i))...
                    + h(j,i) .* (vv_t(1)*e_effx(j,i) + vv_t(2)*e_effy(j,i)) .* (vv_t(1)*spdx(j,i) + vv_t(2)*spdy(j,i))...
                    + h(j,i) .* e_eff(j,i).^(1/3-1) .* (spdxx(j,i).*vv_t(1).^2 + spdxy(j,i).*vv_t(1).*vv_t(2) + spdyx(j,i).*vv_t(1).*vv_t(2) + spdyy(j,i).*vv_t(2).^2));
        angs(j,i) = ang;
        
        bed(j,i) = dr(j,i) + lat(j,i) + lon(j,i);
    end
end

% Masking
dr = dr .* mask;
lat = lat .* mask;
lon = lon .* mask;
bed = bed .* mask;
spd2 = spd2 .* mask;
%%
figure(1)
clf
ax(1) = subplot(121);
p = surf(Xi,Yi,zeros(size(ss)),log10(spd2));
hold on 
utmp = cos(angs);
vtmp = sin(angs);
% quiver(Xi,Yi,utmp,vtmp)
title('Ice Speed')
set(p, 'edgecolor', 'none');
view(2)
axis equal
setFontSize(16);
c = colorbar;
c.Label.String = 'Log_{10} Speed [m/yr]';
ax(2) = subplot(122);
p = surf(Xi,Yi,zeros(size(ss)),b_raw);
hold on
contour(xi,yi,spd2, [10, 10] , 'k:','HandleVisibility','off')
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
title('Bed elevation')
set(p, 'edgecolor', 'none');
colormap(ax(2),icey);
caxis([-2000 500])
view(2)
axis equal
setFontSize(16);
c = colorbar;
c.Label.String = 'Bed Elevation [m]';


figure(2)
clf
ax(1) = subplot(121);
p = surf(Xi,Yi,zeros(size(phi_raw)),phi_raw);
hold on
contour(xi,yi,spd2, [10, 10] , 'k:','HandleVisibility','off')
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
title('Hydro Potential (hydroequalib assumption)')
set(p, 'edgecolor', 'none');
colormap(ax(1),'jet');
caxis([0 700])
view(2)
axis equal
setFontSize(16);
c = colorbar;
c.Label.String = '\Phi [?]';

ax(2) = subplot(122);
[phi_x, phi_y] = gradient(phi,dx,dx);
[b_x, b_y] = gradient(b,dx,dx);
sp = 3;
p = surf(Xi,Yi,zeros(size(phi_raw)),sqrt(phi_x.^2 + phi_y.^2));
hold on
set(p, 'edgecolor', 'none');
quiver(xi(1:sp:end),yi(1:sp:end),-phi_x(1:sp:end,1:sp:end),-phi_y(1:sp:end,1:sp:end));
quiver(xi(1:sp:end),yi(1:sp:end),-b_x(1:sp:end,1:sp:end),-b_y(1:sp:end,1:sp:end));
contour(xi,yi,spd2, [10, 10] , 'k:','HandleVisibility','off')
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
hold off
colorbar
colormap(ax(2),flipud(pink))
caxis([0 .01])
title('Gradients')
legend('|\nabla \Phi|','\nabla \Phi','\nabla Bed')
view(2)
axis equal
setFontSize(16);



%% 
figure(3)
clf
sgtitle('Force Budget (Positive is Along Flow)')
colormap redblue
caxis([-1e5 1e5])
subplot(221)
p = surf(Xi,Yi,zeros(size(ss)),dr);
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
title('Driving Force')
allfig2(p,dr)

subplot(222)
p = surf(Xi,Yi,zeros(size(ss)),lon);
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
title('Longitudinal Stresses')
allfig2(p,lon)

subplot(223)
p = surf(Xi,Yi,zeros(size(ss)),lat);
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
title('Lateral Stresses')
allfig2(p,lat)

subplot(224)
p = surf(Xi,Yi,zeros(size(ss)),bed);
title('Bed Drag')
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
allfig2(p,bed)

setFontSize(16)
%% 
% figure(3)
% clf
% sgtitle('Force Budget Fraction (Positive is Along Flow)')
% colormap redblue
% caxis([-1e5 1e5])
% subplot(221)
% p = surf(Xi,Yi,zeros(size(ss)),dr);
% hold on
% contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
% contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
% contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
% title('Driving Force')
% allfig(p)
% 
% subplot(222)
% p = surf(Xi,Yi,zeros(size(ss)),lon./dr);
% hold on
% contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
% contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
% contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
% title('Longitudinal Stresses')
% allfig(p)
% caxis([-1, 1])
% c = colorbar;
% c.Label.String = '[%]';
% 
% subplot(223)
% p = surf(Xi,Yi,zeros(size(ss)),lat./dr);
% hold on
% contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
% contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
% contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
% title('Lateral Stresses')
% allfig(p)
% caxis([-1, 1])
% c = colorbar;
% c.Label.String = '[%]';
% 
% subplot(224)
% p = surf(Xi,Yi,zeros(size(ss)),bed./dr);
% title('Bed Drag')
% hold on
% contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
% contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
% contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
% allfig(p)
% caxis([-1, 1])
% c = colorbar;
% c.Label.String = '[%]';
% 
% setFontSize(16)
%% 
figure(4)
z = abs(lat)./(abs(lon)+1);
p = surf(Xi,Yi,zeros(size(ss)),z);
title('Ratio of Lat to Lon')
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
allfig2(p,z)
caxis([0, 10])
c = colorbar;
c.Label.String = '[%]';

setFontSize(16)

figure(5)
clf
colormap redblue
p = surf(Xi,Yi,zeros(size(ss)),lat);
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
title('Lateral Stresses')
allfig2(p,lat)

function [] = allfig(p)
set(p, 'edgecolor', 'none');
view(2)
c = colorbar;
c.Label.String = '[Pa]';
% caxis([-1.5e5, 1.5e5])
xlabel('Easting')
ylabel('Northing')
axis equal
end

function [] = allfig2(p,z)
set(p, 'edgecolor', 'none');
view(2)
c = colorbar;
c.Label.String = '[Pa]';
si = std(abs(z),0,'all','omitnan');
me = mean(mean(abs(z),'omitnan'));
caxis([-(me+3*si), (me+3*si)])
% caxis([-1.5e5, 1.5e5])
xlabel('Easting')
ylabel('Northing')
axis equal
end

