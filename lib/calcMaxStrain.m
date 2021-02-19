function [] = calcMaxStrain(north,east,z,dx)
% Function to calculate max strain in a region
%     %% Clean Data (not needed)
%     east(isnan(east)) = nanmean(nanmean(east));
%     north(isnan(north)) = nanmean(nanmean(north));
    %% Smoothing Mask
    if(z~=0)
        north = imgaussfilt(north,z);
        east = imgaussfilt(east,z);
    end
    %% Compute Max Strain
    [n,m] = size(north);
    strain = zeros(n-1,m-1,2,2);
    strain(:,:,1,1) = diff(east(1:end-1,:),1,2);
    strain(:,:,2,2) = diff(north(:,1:end-1));
    strain(:,:,1,2) = .5 .* (diff(north(1:end-1,:),1,2) + diff(east(:,1:end-1)));
    strain(:,:,2,1) =  strain(:,:,1,2);
    maxShear = zeros(n-1,m-1);
    maxShearu = zeros(n-1,m-1);
    maxShearv = zeros(n-1,m-1);
    maxShearu2 = zeros(n-1,m-1);
    maxShearv2 = zeros(n-1,m-1);
     for i = 1:n-1
         for j = 1:m-1
            M = squeeze(strain(i,j,:,:));
            if(sum(sum(isnan(M)))==0)
                [V,D] = eig(M);
                e = diag(D);
                maxShear(i,j) = abs(-e(1)+e(2));
                maxShearu(i,j) = V(1,1);
                maxShearv(i,j) = V(2,1);
                maxShearu2(i,j) = V(1,2);
                maxShearv2(i,j) = V(2,2);
            else
                maxShear(i,j) = NaN;
                maxShearu(i,j) = NaN;
                maxShearv(i,j) = NaN;
                maxShearu2(i,j) = NaN;
                maxShearv2(i,j) = NaN;
            end
         end
     end

    figure
    %contourf(flipud(log(maxShear/45))/log(10),[0,-8:.5:-6])
    XX = (1:(size(north,2)-1))*dx/1000;
    YY = (1:(size(north,1)-1))*dx/1000;
    colormap(flipud(gray))
    surf(XX,YY,log10((maxShear)/dx),'edgecolor','none')
    %caxis([-3 -1.5]);
    view(2)
    xlabel('km');
    ylabel('km');
    title(' Max Shear ')
    c = colorbar;
    c.Location = 'eastoutside';
    c.Label.String = 'log_{10} [s^{-1}]';
    colormap parula
    axis equal
    set(gca,'visible','off','xtick',[])
    
%     figure
%     ss = 1;
%     quiver(XX(1:ss:end),YY(1:ss:end),flipud(maxShearu(1:ss:end,1:ss:end)),flipud(maxShearv(1:ss:end,1:ss:end)));
    %hold on
    %quiver(XX(1:ss:end),YY(1:ss:end),flipud(maxShearu2(1:ss:end,1:ss:end)),flipud(maxShearv2(1:ss:end,1:ss:end)));
    
    %saveas(f,string(year) + 'MaxShearNoAxis.png')  
end