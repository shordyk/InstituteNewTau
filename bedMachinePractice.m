clear
load ll.mat
surfAlong = [0];

% I like to use polarsteric coordinates as they have units of meters, so
% when you draw a box of 1e3 x 1e3 you know its a 1 km square. 

%% This section lets you draw any path, like a radar line, and then see the
% entire square around it, comment out if using section below
[x,y] = ll2ps(Latitude,Longitude); 
if(numel(x)> 1000) %downsample to ensure fast run time.
    x = x(1:round(numel(x)/1000):end);
    y = y(1:round(numel(y)/1000):end);
end
surfAlong = bedmachine_interp('surface',x,y);
bedAlong = bedmachine_interp('bed',x,y);
[xx,yy] = meshgrid(x,y);

%% This section lets you pick your own box to plot, it doesn't draw a radar line
% xmin = -2e5;
% xmax = 2e5;
% ymin = -2e5;
% ymax = 2e5;
% dx   = 10e3;
% [xx,yy] = meshgrid(xmin:dx:xmax,ymin:dx:ymax);

%% Looking up data for the grid
bedHeight = bedmachine_interp('bed',xx,yy);
surfaceHeight = bedmachine_interp('surface',xx,yy);

%% Plotting
figure
imshow("BAS.jpg")
A = imread("BAS.jpg");
title('imshow()')
figure
imagesc(A);
colormap gray
title('imagesc()')

%% 
figure
surf(xx,yy,bedHeight,'edgecolor','none') 
% this plot comes out all black if you don't set edge color to none
hold on
surf(xx,yy,surfaceHeight,'edgecolor','none')

xlabel("Easting [m]")% meters "grid east" of south pole
ylabel("Northing [m]") % meters "grid north" of south pole
zlabel("Elevation ASL [m]") % meters above sea level (ASL)

if(sum(surfAlong)~=0)
plot3(x,y,surfAlong,'r-') %plotting radar line along surface
plot3(x,y,bedAlong,'r-') %plotting radar line along surface
end
colorbar

