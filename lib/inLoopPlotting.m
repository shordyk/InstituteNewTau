figure(1);
    clf
        subplot(131)
            trisurf(t,xy(:,1),xy(:,2),log10(sqrt(u.^2 + v.^2)*3.154E7),...
                'edgecolor','none')
            title('Speed')
            view(2)
            axis equal
            xlabel('X')
            ylabel('Y')
            caxis([1 3.6]);
            colorbar
        subplot(132)
            trisurf(t_c,xy_c(:,1),xy_c(:,2),enhance.^(-nn),...
               'edgecolor','none')
            title("Enhancement " + t_i)
            view(2)
            axis equal
            xlabel('X')
            ylabel('Y')
            colorbar
        subplot(133)
            trisurf(t,xy(:,1),xy(:,2),T_bar(xy(:,1),xy(:,2))-273,...
               'edgecolor','none')
            title("average Temperature ")
            view(2)
            axis equal
            xlabel('X')
            ylabel('Y')
            colorbar
            drawnow
    
    figure(2);
    clf
        subplot(131)
            trisurf(t,xy(:,1),xy(:,2),Br(xy(:,1),xy(:,2)),...
                'edgecolor','none')
            title('Br')
            view(2)
            axis equal
            xlabel('X')
            ylabel('Y')
            colorbar
        subplot(132)
             trisurf(t,xy(:,1),xy(:,2),Pe(xy(:,1),xy(:,2)),...
                'edgecolor','none')
            title('Pe')
            view(2)
            axis equal
            xlabel('X')
            ylabel('Y')
            colorbar
        subplot(133)
             trisurf(t,xy(:,1),xy(:,2),La(xy(:,1),xy(:,2)),...
                'edgecolor','none')
            title('\lambda')
            view(2)
            axis equal
            xlabel('X')
            ylabel('Y')
            colorbar
            drawnow  