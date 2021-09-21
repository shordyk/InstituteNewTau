clear;

[seed_lat, seed_lon] = ps2ll(-9.2212e5,2.5977e5);
% gen_vel_profiles(seed_lat, seed_lon,5);

load institute_antiflow/vel_profile_full.mat
figure
measures('speed','institute ice stream','scalelim',[1 600],'mapwidth',800);
plotm(profile_lat,profile_lon,'k','LineWidth',5)
plotm(seed_lat, seed_lon,'rp')
figure
for j = 1:size(profile_path,2)
    subplot(2,size(profile_path,2),j)
    plot(profile_path(:,j)./1e3,profile_cross(:,j),profile_path(:,j)./1e3,profile_along(:,j),'LineWidth',3)
    ylabel('Speed [m/yr]')
    legend('Across-track','Along-track','Location','NW')
    ylim([-100 250])
    xlim([profile_path(1,j)./1e3 profile_path(end,j)./1e3])
    subplot(2,size(profile_path,2),size(profile_path,2)+j)
    [hice, hbed, hwater] = bedmap2_profile(profile_lat(:,j),profile_lon(:,j));
    ylabel('Elevation [m]')
    xlabel('Along-track Distance [km]')
end
%%
load vel_profiles_paul.mat
figure
mapzoomps('Institute Ice Stream');
measuresps('speed')% [x,y,c] = measures_data('speed',xlim,ylim,'xy');
plotps(profile_lat,profile_lon,'k','LineWidth',5)
hold on
% plotm(seed_lat, seed_lon)
bedmachine('surface','contour',0:50:1000,'k')
caxis([1 1200]);
set(gca,'ColorScale','log')
load measuresColor.mat
colormap(gca, measuresColor)
colorbar
%%
figure
for j = 1:size(profile_path,2)
    subplot(2,size(profile_path,2),j)
    plot(profile_path(:,j)./1e3,profile_cross(:,j),profile_path(:,j)./1e3,profile_along(:,j),'LineWidth',3)
    ylabel('Speed [m/yr]')
    legend('Across-track','Along-track','Location','NW')
    ylim([-100 250])
    xlim([profile_path(1,j)./1e3 profile_path(end,j)./1e3])
    subplot(2,size(profile_path,2),size(profile_path,2)+j)
    [hice, hbed, hwater] = bedmap2_profile(profile_lat(:,j),profile_lon(:,j));
    ylabel('Elevation [m]')
    xlabel('Along-track Distance [km]')
end



function gen_vel_profiles(seed_lat_in, seed_lon_in,n)

    seed_lat = seed_lat_in*ones(1,n)-.2*[0:n-1];
    seed_lon = seed_lon_in*ones(1,n)-.2*[0:n-1];
    for j = 1:length(seed_lat)
        profile_lat_temp = [seed_lat(j)];
        profile_lon_temp = [seed_lon(j)];

        for i = 1:90 %50 previously
            [x1,y1] = ll2ps(profile_lat_temp(1),profile_lon_temp(1));
            vx1 = measures_interp('vx',profile_lat_temp(1),profile_lon_temp(1));
            vy1 = measures_interp('vy',profile_lat_temp(1),profile_lon_temp(1));
            v1 = sqrt(vx1^2 + vy1^2);
            [x2,y2] = ll2ps(profile_lat_temp(end),profile_lon_temp(end));
            vx2 = measures_interp('vx',profile_lat_temp(end),profile_lon_temp(end));
            vy2 = measures_interp('vy',profile_lat_temp(end),profile_lon_temp(end));
            v2 = sqrt(vx2^2 + vy2^2);
            
            stepsize = 800;
            x_temp1 = x1 + stepsize*-(vy1/v1);
            y_temp1 = y1 + stepsize*(vx1/v1);
            [lat_temp1,lon_temp1] = ps2ll(x_temp1,y_temp1);
            x_temp2 = x2 + stepsize*(vy2/v2);
            y_temp2 = y2 + stepsize*-(vx2/v2);
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
    save('vel_profiles_paul.mat','profile_along','profile_cross','profile_lat',...
                            'profile_lon','profile_path')
end