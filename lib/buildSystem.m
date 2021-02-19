% Build sytem for FEM
A = sparse(nE,nN);
B = sparse(nE,nN);
D = sparse(nE,nN);
F = zeros(1,nN);

%Construct discretized gradient operators, and element areas
tau_area = zeros(nE,1);
xy_c = zeros(nE,2);
h_av = zeros(nE,1);
for E = 1:nE  % integration over each element
  nodes = t(E,:);
  xy_c(E,:) = [mean(xy(nodes,1)),mean(xy(nodes,2))];
  xyE = [ones(3,1),xy(nodes,:)];
  Area = abs(det(xyE))/2;
  tau_area(E) = Area;
  h_av(E) = mean(h_re(nodes));
  C = inv(xyE);
  A_E = C(2,:);
  B_E = C(3,:);
  D_E = Area/3*ones(1,3);
  F_E = Area/3*ones(1,3);
  A(E,nodes) = A(E,nodes) + A_E;
  B(E,nodes) = B(E,nodes) + B_E;
  D(E,nodes) = D(E,nodes) + D_E;
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

%% create grid to plot element centers
t_c = delaunay(xy_c(:,1),xy_c(:,2));

%% Driving F estimate for plotting/reference
df = sqrt((A*h_s_init(xy(:,1),xy(:,2))).^2 + ... 
    (B*h_s_init(xy(:,1),xy(:,2))).^2)*rho*g.*h_av;
