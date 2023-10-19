clear;

[seed_lat, seed_lon] = ps2ll(-9.2212e5,2.5977e5);
% changed n to 1 from 3
n = 1; % must be odd
gen_vel_profiles(seed_lat, seed_lon,n);

load vel_profiles_paul_04_13.mat
% figure
% measures('speed','institute ice stream','scalelim',[1 600],'mapwidth',800);
% plotm(profile_lat,profile_lon,'k','LineWidth',3)
% plotm(seed_lat, seed_lon,'rp')
% figure
% for j = 1:size(profile_path,2)
%     subplot(2,size(profile_path,2),j)
%     plot(profile_path(:,j)./1e3,profile_cross(:,j),profile_path(:,j)./1e3,profile_along(:,j),'LineWidth',3)
%     ylabel('Speed [m/yr]')
%     legend('Across-track','Along-track','Location','NW')
%     ylim([-100 250])
%     xlim([profile_path(1,j)./1e3 profile_path(end,j)./1e3])
%     subplot(2,size(profile_path,2),size(profile_path,2)+j)
%     [hice, hbed, hwater] = bedmap2_profile(profile_lat(:,j),profile_lon(:,j));
%     ylabel('Elevation [m]')
%     xlabel('Along-track Distance [km]')
% end
%%
load vel_profiles_paul_04_13.mat
figure(1)
mapzoomps('Institute Ice Stream');
measuresps('speed')% [x,y,c] = measures_data('speed',xlim,ylim,'xy');
plotps(profile_lat,profile_lon,'k','LineWidth',5)
% plotps(profile_lat(90,1),profile_lon(90,1),'ws')
% plotps(profile_lat(90,2),profile_lon(90,2),'ws')
% plotps(profile_lat(90,3),profile_lon(90,3),'ws')
hold on
% plotm(seed_lat, seed_lon)
bedmachine('surface','contour',0:50:1000,'k')
bedmachine('gl','linewidth',3)
xlabel('Contour Plot 1')
caxis([1 1200]);
set(gca,'ColorScale','log')
load measuresColor.mat
colormap(gca, measuresColor)
colorbar
%%
figure(2)
for j = 1:size(profile_path,2)
    subplot(2,size(profile_path,2),j)
    plot(profile_path(:,j)./1e3,profile_cross(:,j),profile_path(:,j)./1e3,profile_along(:,j),'LineWidth',3)
    ylabel('Speed [m/yr]')
    legend('Across-track','Along-track','Location','NW')
    ylim([-100 500])
    xlim([profile_path(1,j)./1e3 profile_path(end,j)./1e3])
    subplot(2,size(profile_path,2),size(profile_path,2)+j)
    [hice, hbed, hwater] = bedmap2_profile(profile_lat(:,j),profile_lon(:,j));
    ylabel('Elevation [m]')
    xlabel('Along-track Distance [km]')
end

%%
figure(3)
mapzoomps('Institute Ice Stream');
xx = xlim;
yy = ylim;
[x,y] = meshgrid(xx(1):1e3:xx(2),yy(1):1e3:yy(2));
err = measures_interp('err',x,y)./measures_interp('speed',x,y);
contourf(x,y,err,[0:.05:1],'edgecolor','none')
view(2)
colorbar
hold on
bedmachine('surface','contour',0:50:1000,'k')
bedmachine('gl','linewidth',3)
title('Percent error in MEaSUREs speed')
caxis([0 1])
plotps(profile_lat,profile_lon,'k','LineWidth',5)

%% Get x and y coordinates of main anti-flow line

% for n=1, lateral and longitduinal coordinates are profile_lat / lon

% Converting to x and y 

[x_line1,y_line1] = ll2ps(profile_lat, profile_lon);



%Making "along path" distance variable along contour line



%disp(along_100)

%newlat = griddedInterpolant(Xi',Yi',lat');
%newlon = griddedInterpolant(Xi', Yi', lon');
%newbase = griddedInterpolant(Xi', Yi', bed'); 

figure(4)
plot(x_line1, y_line1, 'k-')
title('X and Y coord')


%figure
%plot(along_1, newlat(C100X', C100Y'),'LineWidth',1.5)
%title('Force along along 100 spd contour')

%% Create along-track coordinate

along_1 = zeros(size(x_line1));
for i = 2:length(x_line1)
    along_1(i) = along_1(i-1) + sqrt((x_line1(i-1) - x_line1(i))^2 + (y_line1(i-1) - yline_1(i))^2);
end

%% Import basal values


%% Plot basal values along line


%% Calculate longitudinal and lateral forces 


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