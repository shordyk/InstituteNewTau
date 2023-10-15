clear;

[seed_lat, seed_lon] = ps2ll(-9.2212e5,2.5977e5);
% changed n to 1 from 3
n = 1; % must be odd
gen_vel_profiles(seed_lat, seed_lon,n);

load vel_profiles_paul_04_13.mat

%% Import values

%tau = defineTau("ISSM_center_stream");
%tau  = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");
%newtau = squeeze(tau(:,:,21));

%

% Get u,v and xy
load data_strainMesh035ISSM_centeryesAdvectNewBase.mat

%load data_strainMesh035ISSMyesAdvectNewBase.mat

% same xy from define tau, this xy are only used to create tau interpolant
xii   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
yii   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
[Xii, Yii] = ndgrid(xii, yii); 



%% Get x and y coordinates of main anti-flow line
% for n=1, lateral and longitduinal coordinates are profile_lat / lon

% Converting to x and y 
[x_line1,y_line1] = ll2ps(profile_lat, profile_lon);

figure
plot(x_line1, y_line1, 'k-')
title('X and Y coord')
%looks good


%% Calculate temp-dependent viscosity 

%Importing results from mapview run. (u,v. t_bar is mean temp at (x,y)
% Modeling my approach off Paul's forceBalance11

g = 9.81;
dx = 259.6;
dy = 252.23 ;
smth = 4e3;

uinit = u;
vinit = v;

xi = min(xy(:,1))-dx:dx:max(xy(:,1))+dx;
yi = (min(xy(:,2))-dy:dy:max(xy(:,2))+dy)';
%xi = xy(:,1);
%yi = xy(:,2);
[xxx,yyy] = ndgrid(xi,yi);
[Xi,Yi] = meshgrid(xi,yi);

us = scatteredInterpolant(xy(:,1),xy(:,2),u);
vs = scatteredInterpolant(xy(:,1),xy(:,2),v);

enh = scatteredInterpolant(xy_c(:,1), xy_c(:,2), enhance);

v = vs(xxx,yyy);
u = us(xxx,yyy);

spd = sqrt(u.^2 + v.^2);




b_raw =  bedmachine_interp('bed',Xi,Yi);
sf_raw =  bedmachine_interp('surface',Xi,Yi);
b = imgaussfilt(b_raw,2);
sf = imgaussfilt(sf_raw,smth/dx);


% Gradients
h = sf-b;
[ux ,  uy] = gradient(u,dx,dx);
[vx ,  vy] = gradient(v,dx,dx);
[sx ,  sy] = gradient(sf,dx,dx);
[spdx , spdy] = gradient(spd,dx,dx);
[spdxx, spdxy] = gradient(spdx,dx,dx);
[spdyx, spdyy] = gradient(spdy,dx,dx);


% Driving Stress
Tdx = -rho * g * h .* sx.^2;
Tdy = -rho * g * h .* sy.^2;
Td  = sqrt(Tdx.^2 +  Tdy.^2);


% Trying instead to import tau from runs
%tau_c = defineTau("ISSM");
tau_c = defineTau("ISSM_center");
newtau = tau_c(xy(:,1),xy(:,2),uinit,vinit)./norms([uinit,vinit],2,2);
tau_interp = scatteredInterpolant(xy(:,1),xy(:,2),newtau);



%% Calculate viscosity (use . when working with grad? comp?)
% .* for multiplying corresponding components

% Strain rate components

eps_xx = ux;
eps_yy = vy; 
eps_xy = (.5*(uy + vx)).^2;

% Effective strain rate

e_eff = sqrt(.5*(eps_xx.^2 + eps_yy.^2)+ eps_xy);
[e_effx, e_effy] = gradient(e_eff.^(1/3-1),dx,dx);

% Will add "enhance" once this works
taucheck = tau_interp(x_line1, y_line1);
B = 1.6e8;
A = 2.4e-25;
E = enh(xxx,yyy);

%% Calculate longitudinal and lateral forces 

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
        
        lon(j,i) =  2*B*E(j,i)*((vv(1)*sx(j,i) + vv(2)*sy(j,i)) .* e_eff(j,i).^(1/3-1) .* (vv(1)*spdx(j,i) + vv(2)*spdy(j,i))...
                    + h(j,i) .* (vv(1)*e_effx(j,i) + vv(2)*e_effy(j,i)) .* (vv(1)*spdx(j,i) + vv(2)*spdy(j,i))...
                    + h(j,i) .* e_eff(j,i).^(1/3-1) .* (spdxx(j,i).*vv(1).^2 + spdxy(j,i).*vv(1).*vv(2) + spdyx(j,i).*vv(1).*vv(2) + spdyy(j,i).*vv(2).^2))             ;
                
        lat(j,i) =  2*B*E(j,i)*((vv_t(1)*sx(j,i) + vv_t(2)*sy(j,i)) .* e_eff(j,i).^(1/3-1) .* (vv_t(1)*spdx(j,i) + vv_t(2)*spdy(j,i))...
                    + h(j,i) .* (vv_t(1)*e_effx(j,i) + vv_t(2)*e_effy(j,i)) .* (vv_t(1)*spdx(j,i) + vv_t(2)*spdy(j,i))...
                    + h(j,i) .* e_eff(j,i).^(1/3-1) .* (spdxx(j,i).*vv_t(1).^2 + spdxy(j,i).*vv_t(1).*vv_t(2) + spdyx(j,i).*vv_t(1).*vv_t(2) + spdyy(j,i).*vv_t(2).^2));
        angs(j,i) = ang;
        
        bed(j,i) = dr(j,i) + lat(j,i) + lon(j,i);
    end
end

%% Create along-track coordinate

along_1 = zeros(size(x_line1));
for i = 2:length(x_line1)
    along_1(i) = along_1(i-1) + sqrt((x_line1(i-1) - x_line1(i))^2 + (y_line1(i-1) - y_line1(i))^2);
end




%% Lat and Lon along line

lat_interp = griddedInterpolant(xxx,yyy,lat);

lon_interp = griddedInterpolant(xxx,yyy,lon);

%% Mapview Plots - need to edit

figure
clf
sgtitle('Force Budget (Positive is Along Flow)')
colormap redblue
caxis([-1e5 1e5])
subplot(221)
p = surf(xxx,yyy,zeros(size(ss)),dr, "EdgeColor", "none");
hold on
%[C30, H30] = contour(xi,yi,spd, [30, 30] , 'k--','HandleVisibility','off');
%[C100,H] = contour(xi,yi,spd, [100, 100] , 'k-','HandleVisibility','off');
%[C300, H3] = contour(xi,yi,spd, [300, 300] , 'k-','HandleVisibility','off');
%[C3000, H1] = contour(xi,yi,spd, [3000, 3000] , 'k-','HandleVisibility','off');
contour(xi,yi,spd', [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','color',rgb('gray'),'linewidth',2)
title('Driving Force')
colorbar
view(2)


subplot(222)
p = surf(xxx,yyy,zeros(size(ss)),lon, "EdgeColor", "none");
hold on
%contour(xi,yi,spd, [30, 30] , 'k--','HandleVisibility','off')
%contour(xi,yi,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
%contour(xi,yi,spd, [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','color',rgb('gray'),'linewidth',2)
title('Longitudinal Stresses')
colorbar
view(2)

subplot(223)
p = surf(xxx,yyy,zeros(size(ss)),lat, "EdgeColor", "none");
hold on
%contour(xi,yi,spd, [30, 30] , 'k--','HandleVisibility','off')
%contour(xi,yi,spd, [100, 300, 3000] , 'k-','HandleVisibility','off')
%contour(xi,yi,spd, [1000, 1000] , 'k-','LineWidth',2)
bedmachine('gl','color',rgb('gray'),'linewidth',2)
title('Lateral Stresses')
colorbar
view(2)



%% Plot basal values along line

figure
hold on
plot(along_1, lat_interp(x_line1, y_line1),'LineWidth',1.5,'Color', 'g')
plot(along_1, lon_interp(x_line1, y_line1),'LineWidth',1.5,'Color', 'b')
plot(along_1, tau_interp(x_line1, y_line1),'LineWidth',1.5,'Color', 'r' )
legend('Lateral Force', 'Longitudinal Force', 'Basal Drag')
title('Force Balance Along Anti-Flow Line')
hold off


figure
plot(along_1, tau_interp(x_line1, y_line1),'LineWidth',1.5,'Color', 'r' )
title('Basal Drag Along Anti-Flow Line')
%%
function gen_vel_profiles(seed_lat_in, seed_lon_in,n)
    seed_lat = seed_lat_in*ones(1,n)-.4*[-(n-1)/2:(n-1)/2];
    seed_lon = seed_lon_in*ones(1,n)-1.8*[-(n-1)/2:(n-1)/2];
    noise_ratio = .1;
    
    for j = 1:length(seed_lat)
        profile_lat_temp = [seed_lat(j)];
        profile_lon_temp = [seed_lon(j)];
        for i = 1:90 %50 previously
            stepsize = 800;
            err1 = measures_interp('err',profile_lat_temp(1),profile_lon_temp(1));
            speed1 = measures_interp('speed',profile_lat_temp(1),profile_lon_temp(1));
            err2 = measures_interp('err',profile_lat_temp(end),profile_lon_temp(end));
            speed2 = measures_interp('speed',profile_lat_temp(end),profile_lon_temp(end));
            [x1,y1] = ll2ps(profile_lat_temp(1),profile_lon_temp(1));
            if(i < 5 || err1/speed1 < noise_ratio)
                vx1 = measures_interp('vx',profile_lat_temp(1),profile_lon_temp(1));
                vy1 = measures_interp('vy',profile_lat_temp(1),profile_lon_temp(1));
                v1 = sqrt(vx1^2 + vy1^2);
                x_temp1 = x1 + stepsize*-(vy1/v1);
                y_temp1 = y1 + stepsize*(vx1/v1);
            else
%                 disp('high error 1')
                [x_old1, y_old1] = ll2ps(profile_lat_temp(2),profile_lon_temp(2));
                x_temp1 = 2*x1 - x_old1; 
                y_temp1 = 2*y1 - y_old1; 
            end
            [lat_temp1,lon_temp1] = ps2ll(x_temp1,y_temp1);
            [x2,y2] = ll2ps(profile_lat_temp(end),profile_lon_temp(end));
            if(i < 5 || err2/speed2 < noise_ratio)    
                vx2 = measures_interp('vx',profile_lat_temp(end),profile_lon_temp(end));
                vy2 = measures_interp('vy',profile_lat_temp(end),profile_lon_temp(end));
                v2 = sqrt(vx2^2 + vy2^2);
                x_temp2 = x2 + stepsize*(vy2/v2);
                y_temp2 = y2 + stepsize*-(vx2/v2);
            else
%                 disp('high error 2')
                [x_old2, y_old2] = ll2ps(profile_lat_temp(end-1),profile_lon_temp(end-1));
                x_temp2 = 2*x2 - x_old2; 
                y_temp2 = 2*y2 - y_old2;
            end
            
            [lat_temp2,lon_temp2] = ps2ll(x_temp2,y_temp2);
            profile_lat_temp = [lat_temp1,profile_lat_temp,lat_temp2];
            profile_lon_temp = [lon_temp1,profile_lon_temp,lon_temp2];
        end

        profile_lat(:,j) = profile_lat_temp;
        profile_lon(:,j) = profile_lon_temp;
        profile_cross(:,j) = -1*measures_interp('cross',profile_lat_temp,profile_lon_temp);
        profile_along(:,j) = measures_interp('along',profile_lat_temp,profile_lon_temp);
        profile_path(:,j) = pathdist(profile_lat_temp,profile_lon_temp);
    end
    save('vel_profiles_paul_04_13.mat','profile_along','profile_cross','profile_lat',...
                             'profile_lon','profile_path')
end