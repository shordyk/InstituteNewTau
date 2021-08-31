
load('gridInstitute3000.mat')
figure(1)
clf
spd2 = measures_interp('speed',xy(:,1),xy(:,2)); %[m/yr]
trisurf(t,xy(:,1),xy(:,2),zeros(size(spd2)),(spd2),...
       'edgecolor','none')
colorbar
view(2)
axis equal
hold on

xx = -9.4e5:1e4:-8e5;
yy = -sqrt((1.4e5)^2-(xx-(-8.1e5)).^2)+3.85e5;

plot(xx,yy)

figure(2)
clf
trisurf(t,xy(:,1),xy(:,2),heaviside(-(xy(:,1)-(-8.1e5)).^2 - (xy(:,2)-(3.85e5)).^2 + 1.4e5^2),...
       'edgecolor','none','facecolor','interp')
axis equal
view(2)
hold on
colorbar