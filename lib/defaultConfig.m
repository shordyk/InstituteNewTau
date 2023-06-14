function [default] = defaultConfig()
% [default] = defaultConfig() returns default config structure for
% Institute Mapview solver
    default = struct('saveFigs', false,...
                     'saveData', true,...
                     'thermocouple', true,...
                     'plotFigs', true); 
    if(~ismac && ~ispc)
        default.plotFigs = false;
        disp('plotting disabled on linux')
    end
end