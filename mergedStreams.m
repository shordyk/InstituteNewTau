clear;

dx = 1e3;
dy = 1e3;

xmax = 10e5;
ymax = 5e5;
hmax = 3e3;
flow_radius = 5e5;

xi = -xmax/2:dx:xmax/2;
yi = 0:dy:ymax;

[Xi,Yi] = meshgrid(xi,yi);

% Functions
bed_h =@(x,y) subplus((y-50e3)*.002)-500;
surface_h =@(x,y) hmax*(-((y-ymax)/(1.1*ymax)).^2 + (1).^2).^(1/2);
% flow_line = @(y) sqrt((flow_radius)^2-y.^2)-flow_radius;
% flow_line2 = @(y) -sqrt((flow_radius)^2-y.^2)+flow_radius;
flow_line  = @(y)  y/ymax*(xmax/4);
flow_line2 = @(y) -y/ymax*(xmax/4);
tau = @(x,y) max(subplus(100e3 ...
        - 80e3*exp(-(x-flow_line(y)).^2/7e9)- 80e3*exp(-(x-flow_line2(y)).^2/7e9)),...
        subplus((y-ymax)+100e3));

figure(1)
clf
% surf(xi,yi,bed_h(Xi,Yi),'edgecolor', 'none')
% hold on
% surf(xi,yi,surface_h(Xi,Yi),'edgecolor', 'none')
surf(xi,yi,tau(Xi,Yi),'edgecolor', 'none')
hold on
plot(flow_line(yi),yi)
plot(flow_line2(yi),yi)
view(2)