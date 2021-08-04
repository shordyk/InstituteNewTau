clc
clear
%% Initialization
% Physical parameters
H = 1000;
rho = 900;
rho_w = 1000;
g = 9.8;
p = 4/3;
theta = pi/4;

a = 3.7E8/2;
H_w = .8037*H;
alpha = 1E-3;
mu = 2.7E-3;
f = rho*g*sin(alpha);

dx = H*0.057;                     %dx: nominal grid spacing

tau_c =@(y) mu*cos(alpha)*(heaviside(y - H_w).*rho*g.*(H - y) + ...
                           heaviside(H_w - y).*(rho*g*(H - y) - rho_w*g*(H_w/H - y)));


%% Mesh Generation
pv = H*[0,1;2/tan(theta),1;1/tan(theta),0;0,1]; %<- For Schoof benchmark case
[xy,t] = distmesh2d(@dpoly,@huniform,dx,H*[0,0;2/tan(theta),1],pv,pv);
nN = size(xy,1);                     %nN: number of nodes
nE = size(t,1);                      %nE: number of elements
b = unique(boundedges(xy,t));        %b:  boundary node numbers
nB = size(b,1);                      %nB: number of boundary nodes

e = [t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];
e = sort(e,2);
[foo,~,ifoo] = unique(e,'rows');
eB = foo(accumarray(ifoo,1) == 1,:); %eB: boundary edges


%% Build System
A = sparse(nE,nN);
B = sparse(nE,nN);
D = zeros(1,nB);
F = zeros(1,nN);

tau_area = zeros(nE,1);
for E = 1:nE  % integration over each element
  nodes = t(E,:);
  xyE = [ones(3,1),xy(nodes,:)];
  Area = abs(det(xyE))/2;
  tau_area(E) = Area;
  C = inv(xyE);
  A_E = C(2,:);
  B_E = C(3,:);
  F_E = f*Area/3*ones(1,3);
  A(E,nodes) = A(E,nodes) + A_E;
  B(E,nodes) = B(E,nodes) + B_E;
  F(nodes) = F(nodes) + F_E;
end

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
    variables u(nN)    
    obj = 1/p*sum(a*tau_area.*pow_pos(norms([A*u,B*u],2,2),p)) - F*u + sum(b_dx.*tau_c(xy(b,2)).*abs(u(b)));
    subject to
        u >= 0;
    minimize(obj)
cvx_end

%% Visualization
u = u*pi*1e7;
trisurf(t,xy(:,1),xy(:,2),u,u,...
       'edgecolor','none','facecolor','interp');
hold on
%tricontour(t,xy(:,1),xy(:,2),u,...
%           [0,.001,.002,.003,.004,.005,.006,.007,.008,.009,.010,.011,.012,.013,.014,.015,.016]);
plot3(xy(b(u(b) > 1e-6),1),xy(b(u(b) > 1e-6),2),max(u)+0*b(u(b) > 1e-6),'r.','MarkerSize',20)
plot3(xy(b(b_dx==0),1),xy(b(b_dx==0),2),max(u)+0*b(b_dx==0),'b.','MarkerSize',20)
view(2), axis equal, colorbar

figure
plot(u(xy(:,2) > H-dx/2))

surf_vel = u(xy(:,2) > H-dx/2);
%save('data/surf_3_0057.mat','surf_vel')