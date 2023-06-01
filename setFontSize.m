function setFontSize(fontSize)
% Set font size for a figure
set(findall(gcf,'-property','FontSize'),'FontSize',fontSize);
end