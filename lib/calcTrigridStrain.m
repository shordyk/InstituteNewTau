function [int] = calcTriGridStrain(u,v,xy,dx)
    % Calculates strain on triangular nodes, returns on intpolation
    % function
    % elements
    xi = min(xy(:,1))-dx:dx:max(xy(:,1))+dx;
    yi = (min(xy(:,2))-dx:dx:max(xy(:,2))+dx)';
    [Xi,Yi] = ndgrid(xi,yi);
    us = scatteredInterpolant(xy(:,1),xy(:,2),u);
    vs = scatteredInterpolant(xy(:,1),xy(:,2),v);

    %% grid to find strain rates
    ug = us(Xi,Yi); % [m/s]
    vg = vs(Xi,Yi); % [m/s]

    %% Smooth
    ug = imgaussfilt(ug,1);
    vg = imgaussfilt(vg,1);
    
    %% Calc Strain on grid
    [n,m] = size(ug);
    strain = zeros(n-1,m-1,2,2);
    strain(:,:,1,1) = diff(vg(1:end-1,:),1,2)/dx;
    strain(:,:,2,2) = diff(ug(:,1:end-1))/dx;
    strain(:,:,1,2) = .5 .* (diff(ug(1:end-1,:),1,2) + diff(vg(:,1:end-1)))/dx;
    strain(:,:,2,1) =  strain(:,:,1,2);
%%   Total  Effective Strain
    EffStrain = sqrt(.5*(strain(:,:,1,1).^2 + strain(:,:,2,2).^2 +...
            (strain(:,:,1,1) + strain(:,:,2)).^2) + strain(:,:,1,2).^2);
    EffStrain(isnan(EffStrain)) = 1e-14; % this doesn't happen anymore Yay!
%%   Calc Max Shear Strain
%    EffStrain = zeros(n-1,m-1);
%      for i = 1:n-1   
%          for j = 1:m-1
%             M = squeeze(strain(i,j,:,:));
%             if(sum(sum(isnan(M)))==0)
%                 [V,D] = eig(M);
%                 e = diag(D);
%                 EffStrain(i,j) = abs(-e(1)+e(2));
%             else
%                 EffStrain(i,j) = 1e-14;
%             end
%          end
%      end
     %% Interpolate back to triangular grid elements
     int = griddedInterpolant(Xi(2:end,2:end)-dx/2, Yi(2:end,2:end)-dx/2, EffStrain,'linear','nearest'); %Evaluate at middle of cell
end