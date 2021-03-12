function [pv] = makePerimeter(corners,dx)
%  [pv] = makePerimeter(corners,dx) Make file that connects points with linearly spaced side points at
% spacing dx. Useful for mesh outlines of abritraty shapes. Function auto
% wraps the points for you.
    n = size(corners,1);
    corners = [corners; corners(1,:)];
    pv = double.empty(0,2);
    for i = 1:n
        dir = corners(i+1,:)-corners(i,:);
        dir = dir./norm(dir);
        pv_t_1 = (corners(i,1):dx * dir(1):corners(i+1,1))';
        pv_t_2 = (corners(i,2):dx * dir(2):corners(i+1,2))';
        if(dir(1) == 0)
            pv_t_1 = corners(i,1)*ones(size(pv_t_2));
        elseif(dir(2) == 0)
            pv_t_2 = corners(i,2)*ones(size(pv_t_1));
        end
        pv = [pv; [pv_t_1,pv_t_2]];
    end
    
%     figure
%     scatter(corners(:,1),corners(:,2),'ks')
%     hold on 
%     scatter(pv(:,1),pv(:,2),'r*')
end