clear;

[seed_lat, seed_lon] = ps2ll(-9.2212e5,2.5977e5);
% changed n to 1 from 3
n = 1; % must be odd
gen_vel_profiles(seed_lat, seed_lon,n);

load vel_profiles_paul_04_13.mat

% same x and y from define tau 
xi   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","x");
yi   = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_init/strbasemag_AIS_JPL1_ISSM_init.nc","y");
[xx,yy] = ndgrid(xi - 3072000,yi - 3072000);


%% Get x and y coordinates of main anti-flow line
% for n=1, lateral and longitduinal coordinates are profile_lat / lon

% Converting to x and y 
[x_line1,y_line1] = ll2ps(profile_lat, profile_lon);

figure
plot(x_line1, y_line1, 'k-')
title('X and Y coord')
%looks good


%% Create along-track coordinate

along_1 = zeros(size(x_line1));
for i = 2:length(x_line1)
    along_1(i) = along_1(i-1) + sqrt((x_line1(i-1) - x_line1(i))^2 + (y_line1(i-1) - y_line1(i))^2);
end

%% Import basal values

%tau = defineTau("ISSM_center_stream");
tau  = ncread("~/Documents/MATLAB/ISSM/JPL1_ISSM_ctrl/strbasemag_AIS_JPL1_ISSM_ctrl.nc","strbasemag");

newtau = squeeze(tau(:,:,21));

tau_interp = griddedInterpolant(xx,yy,newtau);

%% Plot basal values along line

figure
plot(along_1, tau_interp(x_line1, y_line1),'LineWidth',1.5)
title('Basal Strength Along Anti-Flow Line')

%% Calculate temp-dependent viscosity 

%Importing results from mapview run. (u,v. t_bar is mean temp at (x,y)
% Modeling some of my approach of forceBalance11

load data_strainMesh035ISSMyesAdvectNewBase.mat

g = 9.81;
dx = 2e3;
smth = 4e3;

b_raw =  bedmachine_interp('bed',Xi,Yi);
sf_raw =  bedmachine_interp('surface',Xi,Yi);
b = imgaussfilt(b_raw,2);
sf = imgaussfilt(sf_raw,smth/dx);


[uu, vv] = ndgrid(u,v);

% Gradients
h = sf-b;
[ux ,  uy] = gradient(uu,dx,dx);
[vx ,  vy] = gradient(vv,dx,dx);
[sx ,  sy] = gradient(sf,dx,dx);
[spdx , spdy] = gradient(spd,dx,dx);
[spdxx, spdxy] = gradient(spdx,dx,dx);
[spdyx, spdyy] = gradient(spdy,dx,dx);

% Driving Stress
Tdx = -rho * g * h .* sx.^2;
Tdy = -rho * g * h .* sy.^2;
Td  = sqrt(Tdx.^2 +  Tdy.^2);

% Viscosity


%% Calculate longitudinal and lateral forces 



vel = sqrt(u.^2 + v.^2);


% Importing temperature field for effective viscosity - t_z ? but depth av?



% Angle of incline at each point





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