clc
close all
clear

cd radarData_good
files = dir;
directoryNames = {files.name};
directoryNames = directoryNames(~ismember(directoryNames,{'.','..','.DS_Store'}));
directoryNames = directoryNames(~contains(directoryNames,{'_Layers'}));
cd ..

%20,26 are good examples
for i = [15,20,26]%1:length(directoryNames)
    file = "radarData_good/"+ erase(directoryNames{i}, [".mat"])
    thermalPockets;
    clearvars -except i directoryNames
end