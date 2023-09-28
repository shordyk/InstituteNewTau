% Mapping strain and driving force for Institute Ice Stream
% methods of Van der Veen , also some other plotting to get to know the
% area
clear; close all; 
addpath lib
load Dawn.mat
% load vel_profiles_paul_gl.mat  %antiflow line tracks

icey = cbrewer('div','BrBG',48);
rho = 917;
rho_w = 1000;
g = 9.81;
B = 1.6e8; % A = 2.4e-25 Pa^(-3) s^(-1)
A = 2.4e-25;
overgrab = 0;
xmax = -7e5;
xmin = -12e5;
ymax =  5e5;
ymin =  -1e5;

dx = 2e3;
smth = 4e3;
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

% Internal Deformation expected over a locked bed [m/yr]
% u_int = dr^3 * H * A;

u_int = 2 / 4 *abs(bed).^3 .* h * A * 3.154e7;
%%
figure
clf
ax(1) = subplot(121);
p = surf(Xi,Yi,zeros(size(ss)),log10(spd2));
hold on 
utmp = cos(angs);
vtmp = sin(angs);
% quiver(Xi,Yi,utmp,vtmp)
title('Ice Surface Speed')
set(p, 'edgecolor', 'none');
view(2)
axis equal
setFontSize(16);
c = colorbar;
c.Label.String = 'Log_{10} Speed [m/yr]';
ax(2) = subplot(122);
% p = surf(Xi,Yi,zeros(size(ss)),b_raw);
p = surf(Xi,Yi,b_raw);
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

%%
figure
clf
% ax(1) = subplot(121);
% p = surf(Xi,Yi,zeros(size(phi_raw)),phi_raw);
% hold on
% contour(xi,yi,spd2, [10, 10] , 'k:','HandleVisibility','off')
% contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
% contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
% contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
% title('Hydro Potential (hydroequalib assumption)')
% set(p, 'edgecolor', 'none');
% colormap(ax(1),'jet');
% caxis([0 700])
% view(2)
% axis equal
% setFontSize(16);
% c = colorbar;
% c.Label.String = '\Phi [?]';
% 
% ax(2) = subplot(122);
[phi_x, phi_y] = gradient(phi,dx,dx);
[b_x, b_y] = gradient(b,dx,dx);
sp = 3;
p = surf(Xi,Yi,zeros(size(phi_raw)),sqrt(phi_x.^2 + phi_y.^2));
hold on
set(p, 'edgecolor', 'none');
quiver(xi(1:sp:end),yi(1:sp:end),-phi_x(1:sp:end,1:sp:end),-phi_y(1:sp:end,1:sp:end));
quiver(xi(1:sp:end),yi(1:sp:end),-b_x(1:sp:end,1:sp:end),-b_y(1:sp:end,1:sp:end));
bedmachine('gl','color',rgb('gray'),'linewidth',2)
contour(xi,yi,spd2, [10, 10] , 'k:','HandleVisibility','off')
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')

hold off
colorbar
% colormap(ax(2),flipud(pink))
colormap(flipud(pink))
caxis([0 .01])
title('Gradients of Bed and Hydropotential \Phi at Bed')
legend('|\nabla \Phi|','\nabla \Phi','\nabla Bed','grounding line')
view(2)
axis equal
setFontSize(16);
xlim([-1.2e6 -.7e6]);
ylim([-1e5 5e5]);

%%
% lat and long are 301x251 double
% xi is 251 double
% yi is 301 double so lat[y,x]
%for i = 1:length(xi)
   % for j = 1:length(xi)


%figure



%%
% C is contour matrix
% for first contour coordinates, use C(1,:)
figure
clf
sgtitle('Force Budget (Positive is Along Flow)')
colormap redblue
caxis([-1e5 1e5])
subplot(221)
p = surf(Xi,Yi,zeros(size(ss)),dr);
hold on
[C30, H30] = contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
[C100,H] = contour(xi,yi,spd2, [100, 100] , 'k-','HandleVisibility','off');
[C300, H3] = contour(xi,yi,spd2, [300, 300] , 'k-','HandleVisibility','off');
[C3000, H1] = contour(xi,yi,spd2, [3000, 3000] , 'k-','HandleVisibility','off');
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','color',rgb('gray'),'linewidth',2)
title('Driving Force')
allfig2(p,dr)

subplot(222)
p = surf(Xi,Yi,zeros(size(ss)),lon);
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','color',rgb('gray'),'linewidth',2)
title('Longitudinal Stresses')
allfig2(p,lon)

subplot(223)
p = surf(Xi,Yi,zeros(size(ss)),lat);
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','color',rgb('gray'),'linewidth',2)
title('Lateral Stresses')
allfig2(p,lat)

subplot(224)
p = surf(Xi,Yi,zeros(size(ss)),bed);
title('Bed Drag')
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','color',rgb('gray'),'linewidth',2)
allfig2(p,bed)

setFontSize(16)

%% ---------- Sam Edits ---------

%Contour line x and y values

C100X = C100(1,2:C100(2,1));
C100Y = C100(2,2:C100(2,1));
[X100,Y100] = meshgrid(C100X, C100Y);


%Making "along path" distance variable along contour line

along_100 = zeros(size(C100X));
for i = 2:length(C100X)
    along_100(i) = along_100(i-1) + sqrt((C100X(i-1) - C100X(i))^2 + (C100Y(i-1) - C100Y(i))^2);
end

%disp(along_100)

newlat = griddedInterpolant(Xi',Yi',lat');
newlon = griddedInterpolant(Xi', Yi', lon');
newbase = griddedInterpolant(Xi', Yi', bed'); 

figure
plot(C100X, C100Y, 'k-')
%Looks good. 

figure
plot(along_100, newlat(C100X', C100Y'),'LineWidth',1.5)
title('Force along along 100 spd contour')

hold on 

plot(along_100, newlon(C100X', C100Y'),'LineWidth',1.5)

plot(along_100, newbase(C100X', C100Y'),'LineWidth',1.5)

legend('Lateral Stresses', 'Longitudinal Stresses', 'Bed Drag')
xlabel('Distance along the contour')
ylabel('Force')
hold off
%% 
C300X = C300(1,2:C300(2,1));
C300Y = C300(2,2:C300(2,1));
[X300,Y300] = meshgrid(C300X, C300Y);


%Making "along path" distance variable along contour line

along_300 = zeros(size(C300X));
for i = 2:length(C300X)
    along_300(i) = along_300(i-1) + sqrt((C300X(i-1) - C300X(i))^2 + (C300Y(i-1) - C300Y(i))^2);
end


figure
plot(along_300, newlat(C300X', C300Y'),'LineWidth',1.5)
title('Force along along 300 spd contour')


hold on 

plot(along_300, newlon(C300X', C300Y'),'LineWidth',1.5)

plot(along_300, newbase(C300X', C300Y'),'LineWidth',1.5)

legend('Lateral Stresses', 'Longitudinal Stresses', 'Bed Drag')
xlabel('Distance along the contour')
ylabel('Force')
hold off

%% 1000 Speed Contour

C30X = C30(1,2:C30(2,1));
C30Y = C30(2,2:C30(2,1));
[X30,Y30] = meshgrid(C30X, C30Y);


%Making "along path" distance variable along contour line

along_30 = zeros(size(C30X));
for i = 2:length(C30X)
    along_30(i) = along_30(i-1) + sqrt((C30X(i-1) - C30X(i))^2 + (C30Y(i-1) - C30Y(i))^2);
end


figure
plot(along_30, newlat(C30X', C30Y'),'LineWidth',1.5)
title('Force along along 30 spd contour')


hold on 

plot(along_30, newlon(C30X', C30Y'),'LineWidth',1.5)

plot(along_30, newbase(C30X', C30Y'),'LineWidth',1.5)

xlabel('Distance along the contour')
ylabel('Force')

legend('Lateral Stresses', 'Longitudinal Stresses', 'Bed Drag')
hold off


%%
try
    grd = load("grid/strainMesh035.mat");
    fileLoad = true;
catch
    fileLoad = false;
    warning('grid file not found, skipping plot')
end

if(fileLoad)
    figure
    contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
    hold on
    contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
    contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
    surf(Xi,Yi,-1*ones(size(spd2)),log10(spd2),'edgecolor', 'none');
    caxis([-.9 2.8]);
    plot(grd.xy(grd.nw_bound == 1,1),grd.xy(grd.nw_bound == 1,2),'k','linewidth',2)
    plot(grd.xy(grd.ne_bound == 1,1),grd.xy(grd.ne_bound == 1,2),'k','linewidth',2)
    plot(grd.xy(grd.se_bound == 1,1),grd.xy(grd.se_bound == 1,2),'k','linewidth',2)
    plot(grd.xy(grd.sw_bound == 1,1),grd.xy(grd.sw_bound == 1,2),'k','linewidth',2)
%     contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
%     contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
%     contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
    bedmachine('gl','color',rgb('gray'),'linewidth',2)
    title('Domain')
    view(2)
    axis equal
    setFontSize(16);
    c = colorbar;
    c.Label.String = 'Log_{10} Speed [m/yr]';
end

%%
figure %same as before, but larger plot
clf
colormap redblue
p = surf(Xi,Yi,zeros(size(ss)),lat);
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
title('Lateral Stresses')
bedmachine('gl','color',rgb('grey'),'linewidth',2)
allfig2(p,lat)
setFontSize(20);
%%

figure
clf
p = surf(Xi,Yi,b_raw,(spd2));
hold on
set(p, 'edgecolor', 'none');
p = surf(Xi,Yi,sf_raw,(spd2),'facealpha',0.75);
title('3D plot Elevation of ice surface and bed')
set(p, 'edgecolor', 'none');
colormap(parula);
% caxis([-2000 2000])
% axis equal
setFontSize(16);
c = colorbar;
c.Label.String = 'Ice Surface Speed [m/yr]';
f = gca;
f.ColorScale = 'log';
set(p, 'edgecolor', 'none');
%%

figure
clf

p = surf(Xi,Yi,sf_raw,(spd2),'facealpha',.8);
hold on
contour3(Xi,Yi,sf_raw,[-1000:100:2000],'-','color',rgb('white'))
contour3(Xi,Yi,sf_raw,[-1000:500:2000],'k-')
bedmachine('gl','color',rgb('gray'),'linewidth',2)
title('Elevation of Surface and Contours')
set(p, 'edgecolor', 'none');
colormap(parula);
setFontSize(16);
c = colorbar;
c.Label.String = 'Ice Surface Speed [m/yr]';
f = gca;
f.ColorScale = 'log';
set(p, 'edgecolor', 'none');
view(2)

%%
figure
clf
ax(1) = subplot(121);
p = surf(Xi,Yi,zeros(size(ss)),u_int);
title('Internal Flow (Idealized)')
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
allfig2(p,u_int)
colormap(ax(1),cmocean('tempo'))
caxis([0 100])
c = colorbar;
c.Label.String = '[m/yr]';
setFontSize(16);



ax(2) = subplot(122);
p = surf(Xi,Yi,zeros(size(ss)),(abs(u_int)./spd2*100));
title('Internal Creep Factor')
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
allfig2(p,u_int./spd2*100)
c = colorbar;
colormap(ax(2),cmocean('tempo'))
caxis([0,100]);
c.Label.String = '[%]';
f = gca;
f.ColorScale = 'linear';
setFontSize(16);

%%
figure
clf
p = surf(Xi,Yi,zeros(size(h)),h);
title('Thickness')
hold on
contour(xi,yi,spd2, [10, 10] , 'k:','HandleVisibility','off')
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
allfig2(p,h)
c = colorbar;
c.Label.String = '[m]';
caxis([0 3000])
setFontSize(16);
%%
figure
clf
surf(Xi,Yi,zeros(size(b_raw)),b_raw,'edgecolor','none');
hold on
contour(xi,yi,spd2, [10, 10] , 'k:','HandleVisibility','off')
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
title('Bed elevation')
colormap(flipud(cmocean('turbid')));
caxis([-2000 500])
view(2)
axis equal
setFontSize(16);
c = colorbar;
c.Label.String = 'Bed Elevation [m]';

%%
figure
p = surf(Xi,Yi,zeros(size(ss)),sqrt(sx.^2+sy.^2));
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
[C,hh] = contour(xi,yi,sqrt(sx.^2+sy.^2), [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1], 'r-','HandleVisibility','off');  
clabel(C,hh)
title('Surface Slope, a part of driving stress')
bedmachine('gl','color',rgb('gray'),'linewidth',2)
f = gca;
f.ColorScale = 'log';
view(2)
colorbar
set(p, 'edgecolor', 'none');
caxis([1e-6 1e-2])
setFontSize(16);

%%
figure
surf(Xi,Yi,zeros(size(ss)),measures_interp('err',Xi,Yi)./measures_interp('speed',Xi,Yi),...
    'edgecolor', 'none');
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
view(2)
f = gca;
f.ColorScale = 'linear';
colormap(cmocean('amp'));
caxis([0 1])
colorbar
title('Percent error in speed')
setFontSize(16);

%%
figure 
% from Paul's preprints, was cut from final publication so don't over read into this
p = surf(Xi,Yi,zeros(size(ss)),h./(sqrt(sx.^2+sy.^2)*200e3));
hold on
contour(xi,yi,spd2, [30, 30] , 'k--','HandleVisibility','off')
contour(xi,yi,spd2, [100, 300, 3000] , 'k-','HandleVisibility','off')
contour(xi,yi,spd2, [1000, 1000] , 'k-','LineWidth',2)
[C,hh] = contour(xi,yi,h./(sqrt(sx.^2+sy.^2)*200e3),[.1,.3,1,3,10], 'r-','HandleVisibility','off');  
clabel(C,hh)
title('Perterbation Factor')
bedmachine('gl','k','linewidth',2)
colormap(cmocean('curl'))
f = gca;
f.ColorScale = 'log';
view(2)
colorbar
set(p, 'edgecolor', 'none');
caxis([1e-1 1e1])
setFontSize(16);


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

