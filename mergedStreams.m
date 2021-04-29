clear;

dx = 1e3;
dy = 1e3;
dr = 3e3;

xmax = 100e3;
ymax = 250e3;
hmax = 1.5e3;
flow_radius = 50e3;

% Physical parameters
a = 3.7E8;       %a:     flow parameter pre-factor [Pa s^1/3]
p= 4/3;          %p:     flow parameter power []
g = 10;          %g:     acceleration due to gravity [m/s^2]
rho = 900;       %rho:   density of ice [kg/m^3]
rho_w = 1000;    %rho_w: density of water [kg/m^3]




xi = -xmax/2:dx:xmax/2;
yi = 0:dy:ymax;

[Xi,Yi] = meshgrid(xi,yi);

% Functions
bed_h =@(x,y) subplus((y-50e3)*.002)-500;
surface_h =@(x,y) (hmax-200).*y/ymax + 200;
% flow_line = @(y) sqrt((flow_radius)^2-y.^2)-flow_radius;
% flow_line2 = @(y) -sqrt((flow_radius)^2-y.^2)+flow_radius;
flow_line  = @(y)  y/ymax*(xmax/4);
flow_line2 = @(y) -y/ymax*(xmax/4);
tau =@(x,y) max(subplus(180e3 ...
        - 180e3*exp(-(x-flow_line(y)).^2/1e8)- 180e3*exp(-(x-flow_line2(y)).^2/1e8))+10e3,...
        subplus((y-ymax)*10+200e3));

    
tau_c =@(x,y,u,v) tau(x,y) .*norms([u,v],2,2);

driving_f =@(x,y) (hmax-200)/ymax*rho*g*(surface_h(x,y)-bed_h(x,y));

pv = [-xmax/2,0;xmax/2,0;xmax/2,ymax;-xmax/2,ymax;-xmax/2,0];
[xy,t] = distmesh2d(@dpoly,@huniform,dr,[-xmax/2,0;xmax/2,ymax],pv,pv);    
    
nN = size(xy,1);                     %nN: number of nodes
nE = size(t,1);                      %nE: number of elements
b = unique(boundedges(xy,t));        %b:  boundary node numbers
nB = size(b,1);                      %nB: number of boundary nodes

e = [t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];
e = sort(e,2);
[foo,~,ifoo] = unique(e,'rows');
eB = foo(accumarray(ifoo,1) == 1,:); %eB: boundary edges    

f = driving_f(xy(:,1),xy(:,2));
h = surface_h(xy(:,1),xy(:,2))-bed_h(xy(:,1),xy(:,2));

% calving_force: calving force at ice front
calving_force = .5*(1 - rho/rho_w)*rho*g*min(min(h))^2;
% calving_force = 0;

A = sparse(nE,nN);
B = sparse(nE,nN);
D = zeros(1,nB);
F = zeros(1,nN);

%Construct discretized gradient operators, and element areas
tau_area = zeros(nE,1);
h_av = zeros(nE,1);
for E = 1:nE  % integration over each element
  nodes = t(E,:);
  xyE = [ones(3,1),xy(nodes,:)];
  Area = abs(det(xyE))/2;
  tau_area(E) = Area;
  h_av(E) = mean(h(nodes));
  C = inv(xyE);
  A_E = C(2,:);
  B_E = C(3,:);
  F_E = Area/3*ones(1,3);
  A(E,nodes) = A(E,nodes) + A_E;
  B(E,nodes) = B(E,nodes) + B_E;
  F(nodes) = F(nodes) + F_E;
end

%Construct boundary edge lengths
b_dx = zeros(size(b,1),1);
for N = 1:size(b,1) % integration over each boundary node
    edges = eB(sum(eB == b(N),2) == 1,:);
    for i = 1:size(edges,1)
        edge = edges(i,:);
        b_dx(N) = b_dx(N) + .5*sqrt(sum(diff(xy(edge,:),1).^2,2));
    end
end


%% Solve
cvx_begin
    variables u(nN) v(nN)
    obj = 2*a/p*sum(h_av.*tau_area.*pow_pos(norms([A*u,B*v,1/2*(B*u+A*v)],2,2),p)) + ...
          F*tau_c(xy(:,1),xy(:,2),u,v) + F*(f.*v) + sum(dx*calving_force*v((xy(:,2) < 0E3 + dx/2)));
    subject to
        v <= 0;
        
    minimize(obj)
cvx_end


%%
figure(1)
clf 
trisurf(t,xy(:,1),xy(:,2),v*3.154E7,v*3.154E7,...
       'edgecolor','none','facecolor','interp');
hold on
quiver(xy(:,1),xy(:,2),u*3.154E7,v*3.154E7,2);
view(2)
axis tight
colorbar

figure(2)
clf
% surf(xi,yi,bed_h(Xi,Yi),'edgecolor', 'none')
% hold on
% surf(xi,yi,surface_h(Xi,Yi),'edgecolor', 'none')
surf(xi,yi,tau(Xi,Yi),'edgecolor', 'none')
hold on
plot(flow_line(yi),yi)
plot(flow_line2(yi),yi)
view(2)