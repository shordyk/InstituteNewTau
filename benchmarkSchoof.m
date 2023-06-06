% Benchmark 

clc
clear


%% Initialization
dx = 5E3;                     %dx: nominal grid spacing [m]

% Physical parameters
a = 3.7E8;       %a:     flow parameter pre-factor [Pa s^1/3]
p= 4/3;          %p:     flow parameter power []
g = 10;          %g:     acceleration due to gravity [m/s^2]
rho = 900;       %rho:   density of ice [kg/m^3]
rho_w = 1000;    %rho_w: density of water [kg/m^3]






%% Mesh Generation
pv = [-50E3,-50E3;50E3,-50E3;50E3,250E3;-50E3,250E3;-50E3,-50E3]; %<- For Schoof benchmark case
[xy,t] = distmesh2d(@dpoly,@huniform,dx,[-50E3,-50E3;50E3,250E3],pv,pv);
close
% save('workingGrid');
% load('workingGrid');


nN = size(xy,1);                     %nN: number of nodes
nE = size(t,1);                      %nE: number of elements
b = unique(boundedges(xy,t));        %b:  boundary node numbers
nB = size(b,1);                      %nB: number of boundary nodes

e = [t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];
e = sort(e,2);
[foo,~,ifoo] = unique(e,'rows');
eB = foo(accumarray(ifoo,1) == 1,:); %eB: boundary edges


%tau_c: basal strength profile [Pa]
m = 1;
tau_c =@(x,y,u,v) heaviside(y).*(1-heaviside(y-240e3)).*(subplus(100e3*abs(x/50E3).^m-15e3)+10e3).*norms([u,v],2,2) + ...
                  heaviside(-y).*(2e3).*norms([u,v],2,2)+...
                  heaviside(y-240e3).*100e3.*norms([u,v],2,2);


% Physical fields
%h: ice surface height [m]
h = (450*xy(:,2)/250e3)+250;

%calving_force: calving force at ice front
calving_force = .5*(1 - rho/rho_w)*rho*g*min(min(h))^2;

f = 450/(250e3)*h*g*rho*.1;
% f = 0.8*1e3*ones(size(xy(:,2)));  %f:     gravitational driving force []
% f(xy(:,2) > 0e3) = 7e3;
% f(xy(:,2) > 10e3) = 12e3;
% f(xy(:,2) > 70e3) = 15e3;
% f(xy(:,2) > 190e3) = 12e3;
% f(xy(:,2) > 210e3) = 7e3;
% f(xy(:,2) > 230e3) = 0.8e3;
% f = f*1;
%% Build System
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

for fi = 1
  %f:     gravitational driving force []


  
% u(xy(:,1) > 50e3 -dx/2) ==0;
% u(xy(:,1) < -50e3 +dx/2) ==0;  
%   
%% Solve
cvx_begin
    variables u(nN) v(nN)
    obj = 2*a/p*sum(h_av.*tau_area.*pow_pos(norms([A*u,B*v,1/2*(B*u+A*v)],2,2),p)) + ...
          F*tau_c(xy(:,1),xy(:,2),u,v) + F*(f.*v) + sum(dx*calving_force*v((xy(:,2) < 0E3 + dx/2) & (abs(xy(:,1)) < 50E3 - dx/2)));
    subject to
        v <= 0;
        v(xy(:,2) > 250E3 - dx/2) == 0;
        
    minimize(obj)
cvx_end

%% Visualization
fig3 = figure(3);
clf;
% trisurf(t,xy(:,1),xy(:,2),v*3.154E7,v*3.154E7,...
%        'edgecolor','none','facecolor','interp');
% hold on
% [C,hi] = tricontour(t,xy(:,1),xy(:,2),norms([u,v],2,2)*3.154E7,...
%            [0,10,30,100,300,1000]);
%caxis([0 1]);
colormap gray
hold on
quiver(xy(:,1),xy(:,2),u*3.154E7,v*3.154E7,2);
view(2)
% colorbar
axis tight
drawnow
saveas(fig3,"Df" + fi*1000 +".png");

figure(4), clf;
subplot(131)
trisurf(t,xy(:,1),xy(:,2),tau_c(xy(:,1),xy(:,2),u,v)./norms([u,v],2,2),...
'edgecolor','none','facecolor','interp')%[10e3:10e3:100e3]);
view(2)
colorbar
subplot(132)
trisurf(t,xy(:,1),xy(:,2),f,...
       'edgecolor','none','facecolor','interp')
view(2)
colorbar
subplot(133)
trisurf(t,xy(:,1),xy(:,2),h,...
       'edgecolor','none','facecolor','interp')
view(2)
colorbar
end
   